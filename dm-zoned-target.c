// SPDX-License-Identifier: GPL-2.0-only
/*
 * Copyright (C) 2017 Western Digital Corporation or its affiliates.
 *
 * This file is released under the GPL.
 */

#include "dm-zoned.h"

#include <linux/module.h>


#include <linux/blk_types.h>
#include <linux/sched.h>
#include "linux/kernel.h"
#include <linux/blkdev.h>
#include <linux/blk-cgroup.h>
//#include <linux/rbtree.h>

#define	DM_MSG_PREFIX		"zoned"

#define DMZ_MIN_BIOS		8192

#define DMZ_DEFAULT_BUDGET	(1024 * 1024 * 10)		/* 10G */
#define DMZ_REFILL_TIME		(HZ)			/* 100ms -> 1S */

/*
 * Zone BIO context.
 */
struct dmz_bioctx {
	struct dmz_dev		*dev;
	struct dm_zone		*zone;
	struct bio		*bio;
	refcount_t		ref;
};

/*
 * Cgroup work descriptor.
 */

struct dm_cgrp_data {
	struct blkcg_policy_data pd;
	unsigned int weight;
};

struct dm_cgrp_queued {
        struct delayed_work     refill_work;

	struct list_head	list;

	struct dm_cgrp		*cgrp_data;
 	struct bio_list         queue_list;
	struct dmz_target 	*dmz;

	long long 		check_bw;
	bool			checked;

	unsigned long		check_idle;
};

struct dm_cgrp {
	struct blkg_policy_data pd;

	/* Cgroup identifier*/
	int			css_id;
	char			*css_dir;

	/* Cgroup Weight*/
	unsigned int 		weight;
	
	/* to Cgroup I/O control*/
	long long		budget;
	unsigned int		processed_bw;
	unsigned int		queue_cnt;
	bool			queued;
	
	unsigned long		last_added;
	struct dm_cgrp_queued	*delayed;

	/* Cgroup metadata list*/
	struct list_head	list;

	/* mapping check per Cgroup */
	unsigned int 		*cg_chunks;
	//unsigned int		*cg_zones;
	//unsigned int		*zone_to_chunk;
};

struct dm_cgrp_work {
	struct work_struct	work;
	refcount_t		refcount;

	struct dmz_target	*target;

	struct dm_cgrp		*cgrp_data;
	
	struct bio_list		bio_list;
};



/*
 * Target descriptor.
 */
struct dmz_target {
	struct dm_dev		**ddev;
	unsigned int		nr_ddevs;

	unsigned int		flags;

	/* Zoned block device information */
	struct dmz_dev		*dev;

	/* For metadata handling */
	struct dmz_metadata     *metadata;

	/* For cgroup work */
	struct radix_tree_root	cgrp_rxtree;
	struct workqueue_struct	*cgrp_wq;
	struct mutex		cgrp_lock;

	/* For cloned BIOs to zones */
	struct bio_set		bio_set;

	/* For flush */
	spinlock_t		flush_lock;
	struct bio_list		flush_list;
	struct delayed_work	flush_work;
	struct workqueue_struct *flush_wq;
	
	/* Cgroup Work list head */
	struct dm_cgrp		*cgrp_head;
	struct dm_cgrp_queued	*workon_head;

	/* for budget refill */
	struct mutex		refill_lock;
	struct workqueue_struct	*refill_wq;
	
	long long 		high_bw;
	long long		donate;
	long long 		total_bw;
	long long		top_user_bw;
	unsigned int		high_weight;
	unsigned short		check_point;
};

/*
 * Flush intervals (seconds).
 */
#define DMZ_FLUSH_PERIOD	(10 * HZ)

static struct blkcg_policy blkcg_policy_dmz;

static struct dm_cgrp_data *cpd_to_zcd(struct blkcg_policy_data *cpd)
{
        return cpd ? container_of(cpd, struct dm_cgrp_data, pd) : NULL;
}


static struct dm_cgrp_data *blkcg_to_zcd(struct blkcg *blkcg)
{
        return cpd_to_zcd(blkcg_to_cpd(blkcg, &blkcg_policy_dmz));
}

static struct dm_cgrp *pd_to_zc(struct blkg_policy_data *pd)
{
        return pd ? container_of(pd, struct dm_cgrp, pd) : NULL;
}

static struct dm_cgrp *blkg_to_zc(struct blkcg_gq *blkg)
{
        return pd_to_zc(blkg_to_pd(blkg, &blkcg_policy_dmz));
}


/*
 * Target BIO completion.
 */
static inline void dmz_bio_endio(struct bio *bio, blk_status_t status)
{
	struct dmz_bioctx *bioctx =
		dm_per_bio_data(bio, sizeof(struct dmz_bioctx));

	if (status != BLK_STS_OK && bio->bi_status == BLK_STS_OK)
		bio->bi_status = status;
	if (bioctx->dev && bio->bi_status != BLK_STS_OK)
		bioctx->dev->flags |= DMZ_CHECK_BDEV;

	if (refcount_dec_and_test(&bioctx->ref)) {
		struct dm_zone *zone = bioctx->zone;

		if (zone) {
			if (bio->bi_status != BLK_STS_OK &&
			    bio_op(bio) == REQ_OP_WRITE &&
			    dmz_is_seq(zone))
				set_bit(DMZ_SEQ_WRITE_ERR, &zone->flags);
			dmz_deactivate_zone(zone);
		}
		bio_endio(bio);
	}
}

/*
 * Completion callback for an internally cloned target BIO. This terminates the
 * target BIO when there are no more references to its context.
 */
static void dmz_clone_endio(struct bio *clone)
{
	struct dmz_bioctx *bioctx = clone->bi_private;
	blk_status_t status = clone->bi_status;

	bio_put(clone);
	dmz_bio_endio(bioctx->bio, status);
}

/*
 * Issue a clone of a target BIO. The clone may only partially process the
 * original target BIO.
 */
static int dmz_submit_bio(struct dmz_target *dmz, struct dm_zone *zone,
			  struct bio *bio, sector_t chunk_block,
			  unsigned int nr_blocks)
{
	struct dmz_bioctx *bioctx =
		dm_per_bio_data(bio, sizeof(struct dmz_bioctx));
	struct dmz_dev *dev = zone->dev;
	struct bio *clone;

	if (dev->flags & DMZ_BDEV_DYING)
		return -EIO;

	clone = bio_clone_fast(bio, GFP_NOIO, &dmz->bio_set);
	if (!clone)
		return -ENOMEM;
	bio_set_dev(clone, dev->bdev);
	bioctx->dev = dev;
	clone->bi_iter.bi_sector =
		dmz_start_sect(dmz->metadata, zone) + dmz_blk2sect(chunk_block);
	clone->bi_iter.bi_size = dmz_blk2sect(nr_blocks) << SECTOR_SHIFT;
	clone->bi_end_io = dmz_clone_endio;
	clone->bi_private = bioctx;

	bio_advance(bio, clone->bi_iter.bi_size);

	refcount_inc(&bioctx->ref);
	submit_bio_noacct(clone);
	

	if (bio_op(bio) == REQ_OP_WRITE && dmz_is_seq(zone))
		zone->wp_block += nr_blocks;

	return 0;
}

/*
 * Zero out pages of discarded blocks accessed by a read BIO.
 */
static void dmz_handle_read_zero(struct dmz_target *dmz, struct bio *bio,
				 sector_t chunk_block, unsigned int nr_blocks)
{
	unsigned int size = nr_blocks << DMZ_BLOCK_SHIFT;

	/* Clear nr_blocks */
	swap(bio->bi_iter.bi_size, size);
	zero_fill_bio(bio);
	swap(bio->bi_iter.bi_size, size);

	bio_advance(bio, size);
}

/*
 * Process a read BIO.
 */
static int dmz_handle_read(struct dmz_target *dmz, struct dm_zone *zone,
			   struct bio *bio)
{
	struct dmz_metadata *zmd = dmz->metadata;
	sector_t chunk_block = dmz_chunk_block(zmd, dmz_bio_block(bio));
	unsigned int nr_blocks = dmz_bio_blocks(bio);
	sector_t end_block = chunk_block + nr_blocks;
	struct dm_zone *rzone, *bzone;
	int ret;

	/* Read into unmapped chunks need only zeroing the BIO buffer */
	if (!zone) {
		zero_fill_bio(bio);
		return 0;
	}

	DMDEBUG("(%s): READ chunk %llu -> %s zone %u, block %llu, %u blocks",
		dmz_metadata_label(zmd),
		(unsigned long long)dmz_bio_chunk(zmd, bio),
		(dmz_is_rnd(zone) ? "RND" :
		 (dmz_is_cache(zone) ? "CACHE" : "SEQ")),
		zone->id,
		(unsigned long long)chunk_block, nr_blocks);

	/* Check block validity to determine the read location */
	bzone = zone->bzone;
	while (chunk_block < end_block) {
		nr_blocks = 0;
		if (dmz_is_rnd(zone) || dmz_is_cache(zone) ||
		    chunk_block < zone->wp_block) {
			/* Test block validity in the data zone */
			ret = dmz_block_valid(zmd, zone, chunk_block);
			if (ret < 0)
				return ret;
			if (ret > 0) {
				/* Read data zone blocks */
				nr_blocks = ret;
				rzone = zone;
			}
		}

		/*
		 * No valid blocks found in the data zone.
		 * Check the buffer zone, if there is one.
		 */
		if (!nr_blocks && bzone) {
			ret = dmz_block_valid(zmd, bzone, chunk_block);
			if (ret < 0)
				return ret;
			if (ret > 0) {
				/* Read buffer zone blocks */
				nr_blocks = ret;
				rzone = bzone;
			}
		}

		if (nr_blocks) {
			/* Valid blocks found: read them */
			nr_blocks = min_t(unsigned int, nr_blocks,
					  end_block - chunk_block);
			ret = dmz_submit_bio(dmz, rzone, bio,
					     chunk_block, nr_blocks);
			if (ret)
				return ret;
			chunk_block += nr_blocks;
		} else {
			/* No valid block: zeroout the current BIO block */
			dmz_handle_read_zero(dmz, bio, chunk_block, 1);
			chunk_block++;
		}
	}

	return 0;
}

/*
 * Write blocks directly in a data zone, at the write pointer.
 * If a buffer zone is assigned, invalidate the blocks written
 * in place.
 */
static int dmz_handle_direct_write(struct dmz_target *dmz,
				   struct dm_zone *zone, struct bio *bio,
				   sector_t chunk_block,
				   unsigned int nr_blocks)
{
	struct dmz_metadata *zmd = dmz->metadata;
	struct dm_zone *bzone = zone->bzone;
	int ret;

	if (dmz_is_readonly(zone))
		return -EROFS;

	/* Submit write */
	ret = dmz_submit_bio(dmz, zone, bio, chunk_block, nr_blocks);
	if (ret)
		return ret;

	/*
	 * Validate the blocks in the data zone and invalidate
	 * in the buffer zone, if there is one.
	 */
	ret = dmz_validate_blocks(zmd, zone, chunk_block, nr_blocks);
	if (ret == 0 && bzone)
		ret = dmz_invalidate_blocks(zmd, bzone, chunk_block, nr_blocks);

	return ret;
}

/*
 * Write blocks in the buffer zone of @zone.
 * If no buffer zone is assigned yet, get one.
 * Called with @zone write locked.
 */
static int dmz_handle_buffered_write(struct dmz_target *dmz,
				     struct dm_zone *zone, struct bio *bio,
				     sector_t chunk_block,
				     unsigned int nr_blocks)
{
	struct dmz_metadata *zmd = dmz->metadata;
	struct dm_zone *bzone;
	int ret;

	/* Get the buffer zone. One will be allocated if needed */
	bzone = dmz_get_chunk_buffer(zmd, zone);
	if (IS_ERR(bzone))
		return PTR_ERR(bzone);

	if (dmz_is_readonly(bzone))
		return -EROFS;

	/* Submit write */
	ret = dmz_submit_bio(dmz, bzone, bio, chunk_block, nr_blocks);
	if (ret)
		return ret;

	/*
	 * Validate the blocks in the buffer zone
	 * and invalidate in the data zone.
	 */
	ret = dmz_validate_blocks(zmd, bzone, chunk_block, nr_blocks);
	if (ret == 0 && chunk_block < zone->wp_block)
		ret = dmz_invalidate_blocks(zmd, zone, chunk_block, nr_blocks);

	return ret;
}

/*
 * Process a write BIO.
 */
static int dmz_handle_write(struct dmz_target *dmz, struct dm_zone *zone,
			    struct bio *bio)
{
	struct dmz_metadata *zmd = dmz->metadata;
	sector_t chunk_block = dmz_chunk_block(zmd, dmz_bio_block(bio));
	unsigned int nr_blocks = dmz_bio_blocks(bio);

	if (!zone)
		return -ENOSPC;

	DMDEBUG("(%s): WRITE chunk %llu -> %s zone %u, block %llu, %u blocks",
		dmz_metadata_label(zmd),
		(unsigned long long)dmz_bio_chunk(zmd, bio),
		(dmz_is_rnd(zone) ? "RND" :
		 (dmz_is_cache(zone) ? "CACHE" : "SEQ")),
		zone->id,
		(unsigned long long)chunk_block, nr_blocks);
	if (dmz_is_rnd(zone) || dmz_is_cache(zone) ||
	    chunk_block == zone->wp_block) {
		/*
		 * zone is a random zone or it is a sequential zone
		 * and the BIO is aligned to the zone write pointer:
		 * direct write the zone.
		 */
		return dmz_handle_direct_write(dmz, zone, bio,
					       chunk_block, nr_blocks);
	}

	/*
	 * This is an unaligned write in a sequential zone:
	 * use buffered write.
	 */
	return dmz_handle_buffered_write(dmz, zone, bio, chunk_block, nr_blocks);
}

/*
 * Process a discard BIO.
 */
static int dmz_handle_discard(struct dmz_target *dmz, struct dm_zone *zone,
			      struct bio *bio)
{
	struct dmz_metadata *zmd = dmz->metadata;
	sector_t block = dmz_bio_block(bio);
	unsigned int nr_blocks = dmz_bio_blocks(bio);
	sector_t chunk_block = dmz_chunk_block(zmd, block);
	int ret = 0;

	/* For unmapped chunks, there is nothing to do */
	if (!zone)
		return 0;

	if (dmz_is_readonly(zone))
		return -EROFS;

	DMDEBUG("(%s): DISCARD chunk %llu -> zone %u, block %llu, %u blocks",
		dmz_metadata_label(dmz->metadata),
		(unsigned long long)dmz_bio_chunk(zmd, bio),
		zone->id,
		(unsigned long long)chunk_block, nr_blocks);

	/*
	 * Invalidate blocks in the data zone and its
	 * buffer zone if one is mapped.
	 */
	if (dmz_is_rnd(zone) || dmz_is_cache(zone) ||
	    chunk_block < zone->wp_block)
		ret = dmz_invalidate_blocks(zmd, zone, chunk_block, nr_blocks);
	if (ret == 0 && zone->bzone)
		ret = dmz_invalidate_blocks(zmd, zone->bzone,
					    chunk_block, nr_blocks);
	return ret;
}

static long long dmz_get_bio_size_kb(struct bio *bio)
{
        long long bio_kb;
        bio_kb = dmz_blk2sect(dmz_bio_blocks(bio)) << SECTOR_SHIFT;
        return bio_kb /= 1024;
}

/*
 * Process a BIO.
 */
static void dmz_handle_bio(struct dmz_target *dmz, struct dm_cgrp *cgrp_data, struct bio *bio)
{
	struct dmz_bioctx *bioctx =
		dm_per_bio_data(bio, sizeof(struct dmz_bioctx));
	struct dmz_metadata *zmd = dmz->metadata;
	struct dm_zone *zone;
	int ret;
	long long bio_kb = dmz_get_bio_size_kb(bio);
	unsigned int chunk = dmz_bio_chunk(zmd, bio);
	/* added */
	int css_id = bio->bi_blkg->blkcg->css.id;
	unsigned int mapped_zone = dmz_check_chunk_mapping(zmd, chunk);

	dmz_lock_metadata(zmd);
	
	if (cgrp_data->weight){
		if (cgrp_data->cg_chunks[chunk]) {
			
		       chunk = cgrp_data->cg_chunks[chunk];
		} else {
			if (mapped_zone && bio_op(bio) == 0) goto meta_read;
			printk("%s: Cgroup %d chunk %u has not indirect data", 
					__func__, css_id, chunk);
			while(mapped_zone) {
				chunk++;
				mapped_zone = dmz_check_chunk_mapping(zmd, chunk);
			}
		       	cgrp_data->cg_chunks[dmz_bio_chunk(zmd, bio)] = chunk;
			printk("%s: Cgroup %d chunk %u's indirect data is %u"
					, __func__, css_id, 
					dmz_bio_chunk(zmd, bio), chunk);
		}
	}
meta_read:
	/*
	 * Get the data zone mapping the chunk. There may be no
	 * mapping for read and discard. If a mapping is obtained,
	 + the zone returned will be set to active state.
	 */
	zone = dmz_get_chunk_mapping(zmd, chunk, bio_op(bio));
	
	if (IS_ERR(zone)) {
		ret = PTR_ERR(zone);
		goto out;
	}
	/*
	 * to count how many times bio passed to zone
	 * zone_to_chunk is find the chunk number to chunk id 
	if(zone != NULL){
		if (cgrp_data->zone_to_chunk[zone->id] == 0) 
			cgrp_data->zone_to_chunk[zone->id] = chunk;

		cgrp_data->cg_zones[zone->id] += 1;
	}
	*/

	/* Process the BIO */
	if (zone) {
		dmz_activate_zone(zone);
		bioctx->zone = zone;
		dmz_reclaim_bio_acc(zone->dev->reclaim);
	}
	
	cgrp_data->processed_bw += bio_kb;
			
	switch (bio_op(bio)) {
	case REQ_OP_READ:
		ret = dmz_handle_read(dmz, zone, bio);
		break;
	case REQ_OP_WRITE:
		ret = dmz_handle_write(dmz, zone, bio);
		break;
	case REQ_OP_DISCARD:
	case REQ_OP_WRITE_ZEROES:
		ret = dmz_handle_discard(dmz, zone, bio);
		break;
	default:
		DMERR("(%s): Unsupported BIO operation 0x%x",
		      dmz_metadata_label(dmz->metadata), bio_op(bio));
		ret = -EIO;
	}

	/*
	 * Release the chunk mapping. This will check that the mapping
	 * is still valid, that is, that the zone used still has valid blocks.
	 */
	if (zone)
		dmz_put_chunk_mapping(zmd, zone);
out:
	dmz_bio_endio(bio, errno_to_blk_status(ret));

	dmz_unlock_metadata(zmd);
}

/*
 * Increment a cgroup work reference counter.
 */

static inline void dmz_get_cgrp_work(struct dm_cgrp_work *cgw)
{
	refcount_inc(&cgw->refcount);
}


/*
 * Decrement a cgroup work reference count and
 * free it if it becomes 0.
 */

static void dmz_put_cgrp_work(struct dm_cgrp_work *cgw)
{
	if (refcount_dec_and_test(&cgw->refcount)) {
		WARN_ON(!bio_list_empty(&cgw->bio_list));
		radix_tree_delete(&cgw->target->cgrp_rxtree, cgw->cgrp_data->css_id);
		kfree(cgw);
	}
}


/*
 * Chunk BIO work function.
 */

static void dmz_cgrp_work(struct work_struct *work)
{
	struct dm_cgrp_work *cgw = container_of(work, struct dm_cgrp_work, work);
	struct dmz_target *dmz = cgw->target;
	struct bio *bio;
	
	mutex_lock(&dmz->cgrp_lock);
	while ((bio = bio_list_pop(&cgw->bio_list))) {
		mutex_unlock(&dmz->cgrp_lock);
		dmz_handle_bio(dmz, cgw->cgrp_data, bio);
		mutex_lock(&dmz->cgrp_lock);
		dmz_put_cgrp_work(cgw);
	}
	if(cgw->cgrp_data->queue_cnt == 0) dmz_put_cgrp_work(cgw);
	mutex_unlock(&dmz->cgrp_lock);
}

/*
 * Flush work.
 */
static void dmz_flush_work(struct work_struct *work)
{
	struct dmz_target *dmz = container_of(work, struct dmz_target, flush_work.work);
	struct bio *bio;
	int ret;

	/* Flush dirty metadata blocks */
	ret = dmz_flush_metadata(dmz->metadata);
	if (ret)
		DMDEBUG("(%s): Metadata flush failed, rc=%d",
			dmz_metadata_label(dmz->metadata), ret);

	/* Process queued flush requests */
	while (1) {
		spin_lock(&dmz->flush_lock);
		bio = bio_list_pop(&dmz->flush_list);
		spin_unlock(&dmz->flush_lock);

		if (!bio)
			break;

		dmz_bio_endio(bio, errno_to_blk_status(ret));
	}

	queue_delayed_work(dmz->flush_wq, &dmz->flush_work, DMZ_FLUSH_PERIOD);
}

static struct dm_cgrp *dmz_check_cgrp_list(struct dmz_target *dmz, int css_id)
{
        struct dm_cgrp *cgrp_mbr;
        list_for_each_entry(cgrp_mbr, &dmz->cgrp_head->list, list){
                if (cgrp_mbr->css_id == css_id) return cgrp_mbr;
        }
        return NULL;
}


static void dmz_refill_budget(struct work_struct *work)
{
        struct dm_cgrp_queued *dcq = container_of(work, struct dm_cgrp_queued, refill_work.work);
	struct dm_cgrp_queued *queued_mbr;
        struct dm_cgrp *cgrp_data = dcq->cgrp_data;
        struct dmz_target *dmz = dcq->dmz;
        unsigned int weight = cgrp_data->weight;
	unsigned int workon_weight = 0;
	unsigned int total_weight = 0;
        struct bio *bio;
        long long bio_kb;
	long long donate;
	long long proc_queue = 0;
	long long total_amount = 0;
	bool need_reset = false;
	bool budget_left = false;
	// change working every times 100ms 
	// and make the idle list and check the other cgroup idle and change giving budget
	
	/* 
	 * Checking Top User Weight
	 * And Check Total Weight to Work_on Cgroup
	 */
	list_for_each_entry(queued_mbr,&dmz->workon_head->list,list){
		total_weight += queued_mbr->cgrp_data->weight;
		if (workon_weight < queued_mbr->cgrp_data->weight){
			workon_weight = queued_mbr->cgrp_data->weight;
			dmz->high_weight = workon_weight;
			dmz->top_user_bw = queued_mbr->check_bw;
		}
	}
	dcq->checked = true;
	dcq->check_bw = cgrp_data->processed_bw;
	/*
	 * checking the last period's total bw
	 */
	list_for_each_entry(queued_mbr,&dmz->workon_head->list,list){
		total_amount += queued_mbr->check_bw;
		if (!queued_mbr->checked) break;
		need_reset = true;
	}
	/*
	 * if all Cgroup is checked, flag reset
	 * And save The Total BW and the highest BW
	 * high_bw is the highest Total BW in every period
	 */
	if (need_reset){
		list_for_each_entry(queued_mbr,&dmz->workon_head->list,list){
			queued_mbr->checked = false;
		}
		dmz->total_bw = total_amount;
		if (dmz->high_bw < total_amount)
			dmz->high_bw = total_amount;
	}
	
	cgrp_data->processed_bw = 0;
	if (cgrp_data->budget > 0) budget_left = true;
	
	
	if (dmz->high_weight == weight){
		cgrp_data->budget = LLONG_MAX;
	} else {
		cgrp_data->budget = (dmz->top_user_bw * weight / dmz->high_weight);
	}

	if (cgrp_data->queue_cnt){
        	mutex_lock(&dmz->refill_lock);
	        while ((bio = bio_list_pop(&dcq->queue_list))){
        	        bio_kb = dmz_get_bio_size_kb(bio);
                	cgrp_data->budget -= bio_kb;
			proc_queue += bio_kb;
               	 	mutex_unlock(&dmz->refill_lock);
                	dmz_handle_bio(dmz, cgrp_data,bio);
                	mutex_lock(&dmz->refill_lock);
                	cgrp_data->queue_cnt -= 1;
                	if (cgrp_data->budget <= 0) break;
		}
		mutex_unlock(&dmz->refill_lock);
	}
	if (time_before(dcq->check_idle + (DMZ_REFILL_TIME - (DMZ_REFILL_TIME / 10)),
			       	cgrp_data->last_added)){
		dcq->check_idle = cgrp_data->last_added;
		queue_delayed_work(dmz->refill_wq,
				&cgrp_data->delayed->refill_work,
				DMZ_REFILL_TIME);
	}
	else if (cgrp_data->queue_cnt){
                queue_delayed_work(dmz->refill_wq,
                                &cgrp_data->delayed->refill_work,
                                DMZ_REFILL_TIME);
	}

        else {
		unsigned int zone;
		unsigned int nr_zones = dmz_nr_zones(dmz->metadata);
                cgrp_data->queued = false;
		dcq->checked = false;
		list_del(&cgrp_data->delayed->list);
		donate = (dmz->high_bw * weight / dmz->high_weight);
		total_weight -= weight;
		list_for_each_entry(queued_mbr,&dmz->workon_head->list,list){
			queued_mbr->cgrp_data->budget += 
				(donate * queued_mbr->cgrp_data->weight / total_weight);
		}
		/*
		for(zone = 0 ; zone < nr_zones ; zone++ ){
			if (cgrp_data->cg_zones[zone] > 0 ){
				printk("%s: cgroup %d processed %u times in to zone %u ",
						__func__,
						cgrp_data->css_id,
						cgrp_data->cg_zones[zone],
						zone);
				dmz_report_zones(dmz->metadata, cgrp_data->zone_to_chunk[zone]);
			}
		}
		*/
		
        }
}


/*
 * Get a chunk work and start it to process a new BIO.
 * If the BIO chunk has no work yet, create one.
 */

static int dmz_queue_cgrp_work(struct dmz_target *dmz, struct bio *bio)
{
	int css_id = bio->bi_blkg->blkcg->css.id;

	struct dm_cgrp_work *cgw;
	struct dm_cgrp *cgrp_data;
	struct dm_cgrp_data *zcd;
	unsigned int *cg_chunks;
	//unsigned int *cg_zones;
	//unsigned int *zone_to_chunk;
	unsigned int nr_chunks = dmz_nr_chunks(dmz->metadata);
	unsigned int nr_zones = dmz_nr_zones(dmz->metadata);
	int ret = 0;
	bool ref_added = false;
	long long bio_kb = dmz_get_bio_size_kb(bio);

	mutex_lock(&dmz->cgrp_lock);

	/* Get the BIO chunk work. If one is not active yet, create one */
	cgw = radix_tree_lookup(&dmz->cgrp_rxtree, css_id);
	if (cgw) {
		dmz_get_cgrp_work(cgw);
		ref_added = true;	
	} else {
		/* Create a new cgroup work */
		cgw = kmalloc(sizeof(struct dm_cgrp_work), GFP_NOIO);
		if (unlikely(!cgw)) {
		 	ret = -ENOMEM;
		  	goto out;
		}
		cgrp_data = dmz_check_cgrp_list(dmz, css_id);
                cgw->cgrp_data = cgrp_data;
		INIT_WORK(&cgw->work, dmz_cgrp_work);
		refcount_set(&cgw->refcount, 1);
		cgw->target = dmz;
		/* dm_cgrp init */
		if (cgw->cgrp_data == NULL){
			struct dm_cgrp_queued *delayed;
			delayed = kmalloc(sizeof(struct dm_cgrp_queued), GFP_NOIO);
			cgrp_data = kmalloc(sizeof(struct dm_cgrp),GFP_NOIO);
			cg_chunks = kzalloc(sizeof(unsigned int) * nr_chunks, GFP_KERNEL);
			//cg_zones = kzalloc(sizeof(unsigned int) * nr_zones, GFP_KERNEL);
			//zone_to_chunk = kzalloc(sizeof(unsigned int) * nr_zones, GFP_KERNEL);
			cgw->cgrp_data = cgrp_data;
			cgw->cgrp_data->css_id = css_id;
			cgw->cgrp_data->css_dir = bio->bi_blkg->blkcg->css.cgroup->kn->name;
			cgw->cgrp_data->budget = DMZ_DEFAULT_BUDGET;
			cgw->cgrp_data->queue_cnt = 0;
			cgw->cgrp_data->queued = false;
			cgw->cgrp_data->processed_bw = 0;
			/* added over 3.3 ver */
			cgw->cgrp_data->cg_chunks = cg_chunks;
			//cgw->cgrp_data->cg_zones = cg_zones;
			//cgw->cgrp_data->zone_to_chunk = zone_to_chunk;
			cgw->cgrp_data->delayed = delayed;
			cgw->cgrp_data->delayed->dmz = dmz;
			cgw->cgrp_data->delayed->check_idle = jiffies;
			cgw->cgrp_data->delayed->cgrp_data = cgrp_data;
			cgw->cgrp_data->delayed->check_bw = 0;
			cgw->cgrp_data->delayed->checked = false;
			INIT_DELAYED_WORK(&cgw->cgrp_data->delayed->refill_work, dmz_refill_budget);
			bio_list_init(&cgw->cgrp_data->delayed->queue_list);
			list_add_tail(&cgw->cgrp_data->list, &dmz->cgrp_head->list);
		}
		/* to change weight, so this line */
		zcd = blkcg_to_zcd(bio->bi_blkg->blkcg);
		cgw->cgrp_data->weight = zcd->weight;
		bio_list_init(&cgw->bio_list);
		ret = radix_tree_insert(&dmz->cgrp_rxtree,css_id, cgw);
		if (unlikely(ret)) {
			kfree(cgw);
			goto out;
		}
	}
	if (!cgw->cgrp_data->queued && cgw->cgrp_data->weight) {
		cgw->cgrp_data->queued = true;
		if (dmz->high_weight != 0 && dmz->top_user_bw != 0) 
			cgw->cgrp_data->budget = 
				dmz->top_user_bw * cgw->cgrp_data->weight / dmz->high_weight;
		list_add_tail(&cgw->cgrp_data->delayed->list, &dmz->workon_head->list);
		queue_delayed_work(dmz->refill_wq,
				&cgw->cgrp_data->delayed->refill_work,
				DMZ_REFILL_TIME);
	}
	
	if (cgw->cgrp_data->weight){
		if (cgw->cgrp_data->budget > 0){
			bio_list_add(&cgw->bio_list, bio);
			cgw->cgrp_data->last_added = jiffies;
			cgw->cgrp_data->budget -= bio_kb;
			if (queue_work(dmz->cgrp_wq, &cgw->work))
				dmz_get_cgrp_work(cgw);
		}
		else {
			cgw->cgrp_data->queue_cnt += 1;
			bio_list_add(&cgw->cgrp_data->delayed->queue_list, bio);
			cgw->cgrp_data->last_added = jiffies;
			//if (!ref_added) dmz_get_cgrp_work(cgw);
		}
	} else {	
		bio_list_add(&cgw->bio_list, bio);
		if (queue_work(dmz->cgrp_wq, &cgw->work))
			dmz_get_cgrp_work(cgw);
	}
	

out:
	mutex_unlock(&dmz->cgrp_lock);
	return ret;
}

/*
 * Check if the backing device is being removed. If it's on the way out,
 * start failing I/O. Reclaim and metadata components also call this
 * function to cleanly abort operation in the event of such failure.
 */
bool dmz_bdev_is_dying(struct dmz_dev *dmz_dev)
{
	if (dmz_dev->flags & DMZ_BDEV_DYING)
		return true;

	if (dmz_dev->flags & DMZ_CHECK_BDEV)
		return !dmz_check_bdev(dmz_dev);

	if (blk_queue_dying(bdev_get_queue(dmz_dev->bdev))) {
		dmz_dev_warn(dmz_dev, "Backing device queue dying");
		dmz_dev->flags |= DMZ_BDEV_DYING;
	}

	return dmz_dev->flags & DMZ_BDEV_DYING;
}

/*
 * Check the backing device availability. This detects such events as
 * backing device going offline due to errors, media removals, etc.
 * This check is less efficient than dmz_bdev_is_dying() and should
 * only be performed as a part of error handling.
 */
bool dmz_check_bdev(struct dmz_dev *dmz_dev)
{
	struct gendisk *disk;

	dmz_dev->flags &= ~DMZ_CHECK_BDEV;

	if (dmz_bdev_is_dying(dmz_dev))
		return false;

	disk = dmz_dev->bdev->bd_disk;
	if (disk->fops->check_events &&
	    disk->fops->check_events(disk, 0) & DISK_EVENT_MEDIA_CHANGE) {
		dmz_dev_warn(dmz_dev, "Backing device offline");
		dmz_dev->flags |= DMZ_BDEV_DYING;
	}

	return !(dmz_dev->flags & DMZ_BDEV_DYING);
}

/*
 * Process a new BIO.
 */
static int dmz_map(struct dm_target *ti, struct bio *bio)
{
	struct dmz_target *dmz = ti->private;
	struct dmz_metadata *zmd = dmz->metadata;
	struct dmz_bioctx *bioctx = dm_per_bio_data(bio, sizeof(struct dmz_bioctx));
	sector_t sector = bio->bi_iter.bi_sector;
	unsigned int nr_sectors = bio_sectors(bio);
	sector_t chunk_sector;
	int ret;

	if (dmz_dev_is_dying(zmd))
		return DM_MAPIO_KILL;

	DMDEBUG("(%s): BIO op %d sector %llu + %u => chunk %llu, block %llu, %u blocks",
		dmz_metadata_label(zmd),
		bio_op(bio), (unsigned long long)sector, nr_sectors,
		(unsigned long long)dmz_bio_chunk(zmd, bio),
		(unsigned long long)dmz_chunk_block(zmd, dmz_bio_block(bio)),
		(unsigned int)dmz_bio_blocks(bio));

	if (!nr_sectors && bio_op(bio) != REQ_OP_WRITE)
		return DM_MAPIO_REMAPPED;

	/* The BIO should be block aligned */
	if ((nr_sectors & DMZ_BLOCK_SECTORS_MASK) || (sector & DMZ_BLOCK_SECTORS_MASK))
		return DM_MAPIO_KILL;

	/* Initialize the BIO context */
	bioctx->dev = NULL;
	bioctx->zone = NULL;
	bioctx->bio = bio;
	refcount_set(&bioctx->ref, 1);

	/* Set the BIO pending in the flush list */
	if (!nr_sectors && bio_op(bio) == REQ_OP_WRITE) {
		spin_lock(&dmz->flush_lock);
		bio_list_add(&dmz->flush_list, bio);
		spin_unlock(&dmz->flush_lock);
		mod_delayed_work(dmz->flush_wq, &dmz->flush_work, 0);
		return DM_MAPIO_SUBMITTED;
	}

	/* Split zone BIOs to fit entirely into a zone */
	chunk_sector = sector & (dmz_zone_nr_sectors(zmd) - 1);
	if (chunk_sector + nr_sectors > dmz_zone_nr_sectors(zmd))
		dm_accept_partial_bio(bio, dmz_zone_nr_sectors(zmd) - chunk_sector);

	/* Now ready to handle this BIO */
	ret = dmz_queue_cgrp_work(dmz, bio);
	if (ret) {
		DMDEBUG("(%s): BIO op %d, can't process chunk %llu, err %i",
			dmz_metadata_label(zmd),
			bio_op(bio), (u64)dmz_bio_chunk(zmd, bio),
			ret);
		return DM_MAPIO_REQUEUE;
	}

	return DM_MAPIO_SUBMITTED;
}

/*
 * Get zoned device information.
 */
static int dmz_get_zoned_device(struct dm_target *ti, char *path,
				int idx, int nr_devs)
{
	struct dmz_target *dmz = ti->private;
	struct dm_dev *ddev;
	struct dmz_dev *dev;
	int ret;
	struct block_device *bdev;

	/* Get the target device */
	ret = dm_get_device(ti, path, dm_table_get_mode(ti->table), &ddev);
	if (ret) {
		ti->error = "Get target device failed";
		return ret;
	}

	bdev = ddev->bdev;
	if (bdev_zoned_model(bdev) == BLK_ZONED_NONE) {
		if (nr_devs == 1) {
			ti->error = "Invalid regular device";
			goto err;
		}
		if (idx != 0) {
			ti->error = "First device must be a regular device";
			goto err;
		}
		if (dmz->ddev[0]) {
			ti->error = "Too many regular devices";
			goto err;
		}
		dev = &dmz->dev[idx];
		dev->flags = DMZ_BDEV_REGULAR;
	} else {
		if (dmz->ddev[idx]) {
			ti->error = "Too many zoned devices";
			goto err;
		}
		if (nr_devs > 1 && idx == 0) {
			ti->error = "First device must be a regular device";
			goto err;
		}
		dev = &dmz->dev[idx];
	}
	dev->bdev = bdev;
	dev->dev_idx = idx;
	(void)bdevname(dev->bdev, dev->name);

	dev->capacity = i_size_read(bdev->bd_inode) >> SECTOR_SHIFT;
	if (ti->begin) {
		ti->error = "Partial mapping is not supported";
		goto err;
	}

	dmz->ddev[idx] = ddev;

	return 0;
err:
	dm_put_device(ti, ddev);
	return -EINVAL;
}

/*
 * Cleanup zoned device information.
 */
static void dmz_put_zoned_device(struct dm_target *ti)
{
	struct dmz_target *dmz = ti->private;
	int i;

	for (i = 0; i < dmz->nr_ddevs; i++) {
		if (dmz->ddev[i]) {
			dm_put_device(ti, dmz->ddev[i]);
			dmz->ddev[i] = NULL;
		}
	}
}

static int dmz_fixup_devices(struct dm_target *ti)
{
	struct dmz_target *dmz = ti->private;
	struct dmz_dev *reg_dev, *zoned_dev;
	struct request_queue *q;
	sector_t zone_nr_sectors = 0;
	int i;

	/*
	 * When we have more than on devices, the first one must be a
	 * regular block device and the others zoned block devices.
	 */
	if (dmz->nr_ddevs > 1) {
		reg_dev = &dmz->dev[0];
		if (!(reg_dev->flags & DMZ_BDEV_REGULAR)) {
			ti->error = "Primary disk is not a regular device";
			return -EINVAL;
		}
		for (i = 1; i < dmz->nr_ddevs; i++) {
			zoned_dev = &dmz->dev[i];
			if (zoned_dev->flags & DMZ_BDEV_REGULAR) {
				ti->error = "Secondary disk is not a zoned device";
				return -EINVAL;
			}
			q = bdev_get_queue(zoned_dev->bdev);
			if (zone_nr_sectors &&
			    zone_nr_sectors != blk_queue_zone_sectors(q)) {
				ti->error = "Zone nr sectors mismatch";
				return -EINVAL;
			}
			zone_nr_sectors = blk_queue_zone_sectors(q);
			zoned_dev->zone_nr_sectors = zone_nr_sectors;
			zoned_dev->nr_zones =
				blkdev_nr_zones(zoned_dev->bdev->bd_disk);
		}
	} else {
		reg_dev = NULL;
		zoned_dev = &dmz->dev[0];
		if (zoned_dev->flags & DMZ_BDEV_REGULAR) {
			ti->error = "Disk is not a zoned device";
			return -EINVAL;
		}
		q = bdev_get_queue(zoned_dev->bdev);
		zoned_dev->zone_nr_sectors = blk_queue_zone_sectors(q);
		zoned_dev->nr_zones = blkdev_nr_zones(zoned_dev->bdev->bd_disk);
	}

	if (reg_dev) {
		sector_t zone_offset;

		reg_dev->zone_nr_sectors = zone_nr_sectors;
		reg_dev->nr_zones =
			DIV_ROUND_UP_SECTOR_T(reg_dev->capacity,
					      reg_dev->zone_nr_sectors);
		reg_dev->zone_offset = 0;
		zone_offset = reg_dev->nr_zones;
		for (i = 1; i < dmz->nr_ddevs; i++) {
			dmz->dev[i].zone_offset = zone_offset;
			zone_offset += dmz->dev[i].nr_zones;
		}
	}
	return 0;
}

/*
 * Setup target.
 */
static int dmz_ctr(struct dm_target *ti, unsigned int argc, char **argv)
{
	struct dmz_target *dmz;
	struct dm_cgrp *cgrp_head;
	struct dm_cgrp_queued *workon_head;
	int ret, i;

	/* Check arguments */
	if (argc < 1) {
		ti->error = "Invalid argument count";
		return -EINVAL;
	}

	/* Allocate and initialize the target descriptor */
	dmz = kzalloc(sizeof(struct dmz_target), GFP_KERNEL);
	if (!dmz) {
		ti->error = "Unable to allocate the zoned target descriptor";
		return -ENOMEM;
	}
	dmz->dev = kcalloc(argc, sizeof(struct dmz_dev), GFP_KERNEL);
	if (!dmz->dev) {
		ti->error = "Unable to allocate the zoned device descriptors";
		kfree(dmz);
		return -ENOMEM;
	}
	dmz->ddev = kcalloc(argc, sizeof(struct dm_dev *), GFP_KERNEL);
	if (!dmz->ddev) {
		ti->error = "Unable to allocate the dm device descriptors";
		ret = -ENOMEM;
		goto err;
	}
	dmz->nr_ddevs = argc;

	ti->private = dmz;

	/* Get the target zoned block device */
	for (i = 0; i < argc; i++) {
		ret = dmz_get_zoned_device(ti, argv[i], i, argc);
		if (ret)
			goto err_dev;
	}
	ret = dmz_fixup_devices(ti);
	if (ret)
		goto err_dev;

	/* Initialize metadata */
	ret = dmz_ctr_metadata(dmz->dev, argc, &dmz->metadata,
			       dm_table_device_name(ti->table));
	if (ret) {
		ti->error = "Metadata initialization failed";
		goto err_dev;
	}

	/* Set target (no write same support) */
	ti->max_io_len = dmz_zone_nr_sectors(dmz->metadata);
	ti->num_flush_bios = 1;
	ti->num_discard_bios = 1;
	ti->num_write_zeroes_bios = 1;
	ti->per_io_data_size = sizeof(struct dmz_bioctx);
	ti->flush_supported = true;
	ti->discards_supported = true;

	/* The exposed capacity is the number of chunks that can be mapped */
	ti->len = (sector_t)dmz_nr_chunks(dmz->metadata) <<
		dmz_zone_nr_sectors_shift(dmz->metadata);

	/* Zone BIO */
	ret = bioset_init(&dmz->bio_set, DMZ_MIN_BIOS, 0, 0);
	if (ret) {
		ti->error = "Create BIO set failed";
		goto err_meta;
	}

	/* Chunk BIO work */
	mutex_init(&dmz->cgrp_lock);
	INIT_RADIX_TREE(&dmz->cgrp_rxtree, GFP_NOIO);
	dmz->cgrp_wq = alloc_workqueue("dmz_cgwq_%s",
					WQ_MEM_RECLAIM | WQ_UNBOUND, 0,
					dmz_metadata_label(dmz->metadata));

	if (!dmz->cgrp_wq) {
		ti->error = "Create cgroup workqueue failed";
		ret = -ENOMEM;
                goto err_bio;
        }

	mutex_init(&dmz->refill_lock);
	dmz->refill_wq = alloc_ordered_workqueue("dmz_rfwq_%s", WQ_MEM_RECLAIM,
						dmz_metadata_label(dmz->metadata));

	/* Flush work */
	spin_lock_init(&dmz->flush_lock);
	bio_list_init(&dmz->flush_list);
	INIT_DELAYED_WORK(&dmz->flush_work, dmz_flush_work);
	dmz->flush_wq = alloc_ordered_workqueue("dmz_fwq_%s", WQ_MEM_RECLAIM,
						dmz_metadata_label(dmz->metadata));
	if (!dmz->flush_wq) {
		ti->error = "Create flush workqueue failed";
		ret = -ENOMEM;
		goto err_cwq;
	}
	mod_delayed_work(dmz->flush_wq, &dmz->flush_work, DMZ_FLUSH_PERIOD);

	/* Initialize reclaim */
	for (i = 0; i < dmz->nr_ddevs; i++) {
		ret = dmz_ctr_reclaim(dmz->metadata, &dmz->dev[i].reclaim, i);
		if (ret) {
			ti->error = "Zone reclaim initialization failed";
			goto err_fwq;
		}
	}

	/*
	 * dm_cgrp_head init
	 */
	cgrp_head = kmalloc(sizeof(struct dm_cgrp), GFP_NOIO);
	dmz->cgrp_head = cgrp_head;
	dmz->cgrp_head->css_id = 0;
	dmz->cgrp_head->css_dir = "HEAD";
	dmz->cgrp_head->weight = 0;
	INIT_LIST_HEAD(&dmz->cgrp_head->list);
	
	/*
	 *
	 */
	workon_head = kmalloc(sizeof(struct dm_cgrp_queued), GFP_NOIO);
	dmz->workon_head = workon_head;
	INIT_LIST_HEAD(&dmz->workon_head->list);

	dmz->high_weight = 0;
	dmz->high_bw = 0;
	dmz->donate = 0;
	dmz->top_user_bw = DMZ_DEFAULT_BUDGET;
	dmz->total_bw = DMZ_DEFAULT_BUDGET;
	dmz->check_point = 0;

	DMINFO("(%s): Target device: %llu 512-byte logical sectors (%llu blocks)",
	       dmz_metadata_label(dmz->metadata),
	       (unsigned long long)ti->len,
	       (unsigned long long)dmz_sect2blk(ti->len));

	return 0;
err_fwq:
	destroy_workqueue(dmz->flush_wq);
err_cwq:
	destroy_workqueue(dmz->cgrp_wq);
err_bio:
	mutex_destroy(&dmz->cgrp_lock);
	bioset_exit(&dmz->bio_set);
err_meta:
	dmz_dtr_metadata(dmz->metadata);
err_dev:
	dmz_put_zoned_device(ti);
err:
	kfree(dmz->dev);
	kfree(dmz);

	return ret;
}

/*
 * Cleanup target.
 */
static void dmz_dtr(struct dm_target *ti)
{
	struct dmz_target *dmz = ti->private;
	struct dm_cgrp *cgrp_mbr;
	int i;

	flush_workqueue(dmz->cgrp_wq);
	destroy_workqueue(dmz->cgrp_wq);
	
	list_for_each_entry(cgrp_mbr, &dmz->cgrp_head->list, list) kfree(cgrp_mbr);
	kfree(dmz->cgrp_head);

	for (i = 0; i < dmz->nr_ddevs; i++)
		dmz_dtr_reclaim(dmz->dev[i].reclaim);

	cancel_delayed_work_sync(&dmz->flush_work);
	destroy_workqueue(dmz->flush_wq);

	(void) dmz_flush_metadata(dmz->metadata);

	dmz_dtr_metadata(dmz->metadata);

	bioset_exit(&dmz->bio_set);

	dmz_put_zoned_device(ti);

	mutex_destroy(&dmz->cgrp_lock);

	kfree(dmz->dev);
	kfree(dmz);
}

/*
 * Setup target request queue limits.
 */
static void dmz_io_hints(struct dm_target *ti, struct queue_limits *limits)
{
	struct dmz_target *dmz = ti->private;
	unsigned int chunk_sectors = dmz_zone_nr_sectors(dmz->metadata);

	limits->logical_block_size = DMZ_BLOCK_SIZE;
	limits->physical_block_size = DMZ_BLOCK_SIZE;

	blk_limits_io_min(limits, DMZ_BLOCK_SIZE);
	blk_limits_io_opt(limits, DMZ_BLOCK_SIZE);

	limits->discard_alignment = DMZ_BLOCK_SIZE;
	limits->discard_granularity = DMZ_BLOCK_SIZE;
	limits->max_discard_sectors = chunk_sectors;
	limits->max_hw_discard_sectors = chunk_sectors;
	limits->max_write_zeroes_sectors = chunk_sectors;

	/* FS hint to try to align to the device zone size */
	limits->chunk_sectors = chunk_sectors;
	limits->max_sectors = chunk_sectors;

	/* We are exposing a drive-managed zoned block device */
	limits->zoned = BLK_ZONED_NONE;
}

/*
 * Pass on ioctl to the backend device.
 */
static int dmz_prepare_ioctl(struct dm_target *ti, struct block_device **bdev)
{
	struct dmz_target *dmz = ti->private;
	struct dmz_dev *dev = &dmz->dev[0];

	if (!dmz_check_bdev(dev))
		return -EIO;

	*bdev = dev->bdev;

	return 0;
}

/*
 * Stop works on suspend.
 */
static void dmz_suspend(struct dm_target *ti)
{
	struct dmz_target *dmz = ti->private;
	int i;

	flush_workqueue(dmz->cgrp_wq);
	for (i = 0; i < dmz->nr_ddevs; i++)
		dmz_suspend_reclaim(dmz->dev[i].reclaim);
	cancel_delayed_work_sync(&dmz->flush_work);
}

/*
 * Restart works on resume or if suspend failed.
 */
static void dmz_resume(struct dm_target *ti)
{
	struct dmz_target *dmz = ti->private;
	int i;

	queue_delayed_work(dmz->flush_wq, &dmz->flush_work, DMZ_FLUSH_PERIOD);
	for (i = 0; i < dmz->nr_ddevs; i++)
		dmz_resume_reclaim(dmz->dev[i].reclaim);
}

static int dmz_iterate_devices(struct dm_target *ti,
			       iterate_devices_callout_fn fn, void *data)
{
	struct dmz_target *dmz = ti->private;
	unsigned int zone_nr_sectors = dmz_zone_nr_sectors(dmz->metadata);
	sector_t capacity;
	int i, r;

	for (i = 0; i < dmz->nr_ddevs; i++) {
		capacity = dmz->dev[i].capacity & ~(zone_nr_sectors - 1);
		r = fn(ti, dmz->ddev[i], 0, capacity, data);
		if (r)
			break;
	}
	return r;
}

static void dmz_status(struct dm_target *ti, status_type_t type,
		       unsigned int status_flags, char *result,
		       unsigned int maxlen)
{
	struct dmz_target *dmz = ti->private;
	ssize_t sz = 0;
	char buf[BDEVNAME_SIZE];
	struct dmz_dev *dev;
	int i;

	switch (type) {
	case STATUSTYPE_INFO:
		DMEMIT("%u zones %u/%u cache",
		       dmz_nr_zones(dmz->metadata),
		       dmz_nr_unmap_cache_zones(dmz->metadata),
		       dmz_nr_cache_zones(dmz->metadata));
		for (i = 0; i < dmz->nr_ddevs; i++) {
			/*
			 * For a multi-device setup the first device
			 * contains only cache zones.
			 */
			if ((i == 0) &&
			    (dmz_nr_cache_zones(dmz->metadata) > 0))
				continue;
			DMEMIT(" %u/%u random %u/%u sequential",
			       dmz_nr_unmap_rnd_zones(dmz->metadata, i),
			       dmz_nr_rnd_zones(dmz->metadata, i),
			       dmz_nr_unmap_seq_zones(dmz->metadata, i),
			       dmz_nr_seq_zones(dmz->metadata, i));
		}
		break;
	case STATUSTYPE_TABLE:
		dev = &dmz->dev[0];
		format_dev_t(buf, dev->bdev->bd_dev);
		DMEMIT("%s", buf);
		for (i = 1; i < dmz->nr_ddevs; i++) {
			dev = &dmz->dev[i];
			format_dev_t(buf, dev->bdev->bd_dev);
			DMEMIT(" %s", buf);
		}
		break;
	}
	return;
}

static int dmz_message(struct dm_target *ti, unsigned int argc, char **argv,
		       char *result, unsigned int maxlen)
{
	struct dmz_target *dmz = ti->private;
	int r = -EINVAL;

	if (!strcasecmp(argv[0], "reclaim")) {
		int i;

		for (i = 0; i < dmz->nr_ddevs; i++)
			dmz_schedule_reclaim(dmz->dev[i].reclaim);
		r = 0;
	} else
		DMERR("unrecognized message %s", argv[0]);
	return r;
}

static struct target_type dmz_type = {
	.name		 = "zoned",
	.version	 = {2, 0, 0},
	.features	 = DM_TARGET_SINGLETON | DM_TARGET_MIXED_ZONED_MODEL,
	.module		 = THIS_MODULE,
	.ctr		 = dmz_ctr,
	.dtr		 = dmz_dtr,
	.map		 = dmz_map,
	.io_hints	 = dmz_io_hints,
	.prepare_ioctl	 = dmz_prepare_ioctl,
	.postsuspend	 = dmz_suspend,
	.resume		 = dmz_resume,
	.iterate_devices = dmz_iterate_devices,
	.status		 = dmz_status,
	.message	 = dmz_message,
};

static int __init dmz_init(void)
{
	blkcg_policy_register(&blkcg_policy_dmz);
	return dm_register_target(&dmz_type);
}

static void __exit dmz_exit(void)
{
	blkcg_policy_unregister(&blkcg_policy_dmz);
	dm_unregister_target(&dmz_type);
}

static int dmz_io_set_weight_legacy(struct cgroup_subsys_state *css,
		struct cftype *cftype,
		u64 val)
{
	struct blkcg *blkcg = css_to_blkcg(css);
	struct dm_cgrp_data *zcd = blkcg_to_zcd(blkcg);
	struct blkcg_gq *blkg;
	int ret = -ERANGE;

	if (val < DMZ_MIN_WEIGHT || val > DMZ_MAX_WEIGHT)
		return ret;
	ret = 0;

	spin_lock_irq(&blkcg->lock);
	zcd->weight = (unsigned int)val;
	hlist_for_each_entry(blkg, &blkcg->blkg_list, blkcg_node) {
		struct dm_cgrp *cgrp_data = blkg_to_zc(blkg);

		if (cgrp_data)
			cgrp_data->weight = zcd->weight;
	}
	spin_unlock_irq(&blkcg->lock);

	return ret;
}

static ssize_t dmz_io_set_weight(struct kernfs_open_file *of,
                char *buf, size_t nbytes,
                loff_t off)
{
        u64 weight;
        int ret = kstrtoull(strim(buf), 0, &weight);

        if (ret)
                return ret;

        ret = dmz_io_set_weight_legacy(of_css(of), NULL, weight);
        return ret ?: nbytes;
}

static int dmz_io_show_weight(struct seq_file *sf, void *v)
{
	struct blkcg *blkcg = css_to_blkcg(seq_css(sf));
	struct dm_cgrp_data *zcd = blkcg_to_zcd(blkcg);
	unsigned int val = 0;

	if (blkcg)
		val = zcd->weight;
	seq_printf(sf, "%u", val);

	return 0;
}

static struct cftype dmz_blkg_files[] = {
	{
		.name = "dmz.weight",
		.flags = CFTYPE_NOT_ON_ROOT,
		.seq_show = dmz_io_show_weight,
		.write = dmz_io_set_weight,
	},
	{ } /* terminate */ 
};

static struct cftype dmz_blkcg_legacy_files[] = {
	{
		.name = "dmz.weight",
		.flags = CFTYPE_NOT_ON_ROOT,
		.seq_show = dmz_io_show_weight,
		.write_u64 = dmz_io_set_weight_legacy,
	},
	{ } /* terminate */
};

static struct blkcg_policy_data *dmz_cpd_alloc(gfp_t gfp)
{
	struct dm_cgrp_data *zcd;

	zcd = kzalloc(sizeof(*zcd), gfp);
	if (!zcd)
		return NULL;
	return &zcd->pd;
}

static void dmz_cpd_init(struct blkcg_policy_data *cpd)
{
	struct dm_cgrp_data *zcd = cpd_to_zcd(cpd);
	zcd->weight = cgroup_subsys_on_dfl(io_cgrp_subsys) ?
		CGROUP_WEIGHT_DFL : DMZ_WEIGHT_LEGACY_DFL;
}

static void dmz_cpd_free(struct blkcg_policy_data *cpd)
{
	kfree(cpd_to_zcd(cpd));
}

static struct blkg_policy_data *dmz_pd_alloc(gfp_t gfp, struct request_queue *q,
						struct blkcg *blkcg)
{
	struct dm_cgrp *zc;

	zc = kzalloc_node(sizeof(*zc), gfp, q->node);
	if (!zc)
		return NULL;
	return &zc->pd;
}

static void dmz_pd_init(struct blkg_policy_data *pd)
{
	struct blkcg_gq *blkg = pd_to_blkg(pd);
	struct dm_cgrp *cgrp_data = blkg_to_zc(blkg);
	struct dm_cgrp_data *zcd = blkcg_to_zcd(blkg->blkcg);

	cgrp_data->weight = zcd->weight;
}

static void dmz_pd_free(struct blkg_policy_data *pd)
{
	struct dm_cgrp *cgrp_data = pd_to_zc(pd);
	kfree(cgrp_data);
}

static struct blkcg_policy blkcg_policy_dmz = {
	.dfl_cftypes		= dmz_blkg_files,
	.legacy_cftypes		= dmz_blkcg_legacy_files,

	.cpd_alloc_fn		= dmz_cpd_alloc,
	.cpd_init_fn		= dmz_cpd_init,
	.cpd_free_fn		= dmz_cpd_free,

	.pd_alloc_fn		= dmz_pd_alloc,
	.pd_init_fn		= dmz_pd_init,
	.pd_free_fn		= dmz_pd_free,
};

module_init(dmz_init);
module_exit(dmz_exit);


MODULE_DESCRIPTION(DM_NAME " target for zoned block devices");
MODULE_AUTHOR("Damien Le Moal <damien.lemoal@wdc.com>");
MODULE_LICENSE("GPL");
