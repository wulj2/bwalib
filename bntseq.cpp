/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <unistd.h>
#include <errno.h>
#include "bntseq.h"
#include "utils.h"
#include "khash.h"
KHASH_MAP_INIT_STR(str, int)

unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void bns_dump(const bntseq_t *bns, const char *prefix)
{
	char str[1024];
	FILE *fp;
	int i;
	{ // dump .ann
		strcpy(str, prefix); strcat(str, ".ann");
		fp = xopen(str, "w");
		err_fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->seed);
		for (i = 0; i != bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			err_fprintf(fp, "%d %s", p->gi, p->name);
			if (p->anno[0]) err_fprintf(fp, " %s\n", p->anno);
			else err_fprintf(fp, "\n");
			err_fprintf(fp, "%lld %d %d\n", (long long)p->offset, p->len, p->n_ambs);
		}
		err_fflush(fp);
		err_fclose(fp);
	}
	{ // dump .amb
		strcpy(str, prefix); strcat(str, ".amb");
		fp = xopen(str, "w");
		err_fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->n_holes);
		for (i = 0; i != bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			err_fprintf(fp, "%lld %d %c\n", (long long)p->offset, p->len, p->amb);
		}
		err_fflush(fp);
		err_fclose(fp);
	}
}

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[8192];
	FILE *fp;
	const char *fname;
	bntseq_t *bns;
	long long xx;
	int i;
	int scanres;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = xopen(fname = ann_filename, "r");
		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		if (scanres != 3) goto badread;
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			scanres = fscanf(fp, "%u%s", &p->gi, str);
			if (scanres != 2) goto badread;
			p->name = strdup(str);
			// read fasta comments 
			while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			while (c != '\n' && c != EOF) c = fgetc(fp);
			if (c == EOF) {
				scanres = EOF;
				goto badread;
			}
			*q = 0;
			if (q - str > 1 && strcmp(str, " (null)") != 0) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			if (scanres != 3) goto badread;
			p->offset = xx;
		}
		err_fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = xopen(fname = amb_filename, "r");
		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		if (scanres != 3) goto badread;
		l_pac = xx;
		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = bns->n_holes? (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t)) : 0;
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			scanres = fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			if (scanres != 3) goto badread;
			p->offset = xx;
			p->amb = str[0];
		}
		err_fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = xopen(pac_filename, "rb");
	}
	return bns;

 badread:
	if (EOF == scanres) {
		err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
	}
	err_fatal(__func__, "Parse error reading %s\n", fname);
}

bntseq_t *bns_restore(const char *prefix)
{  
	char ann_filename[1024], amb_filename[1024], pac_filename[1024], alt_filename[1024];
	FILE *fp;
	bntseq_t *bns;
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	bns = bns_restore_core(ann_filename, amb_filename, pac_filename);
	if (bns == 0) return 0;
	if ((fp = fopen(strcat(strcpy(alt_filename, prefix), ".alt"), "r")) != 0) { // read .alt file if present
		char str[1024];
		khash_t(str) *h;
		int c, i, absent;
		khint_t k;
		h = kh_init(str);
		for (i = 0; i < bns->n_seqs; ++i) {
			k = kh_put(str, h, bns->anns[i].name, &absent);
			kh_val(h, k) = i;
		}
		i = 0;
		while ((c = fgetc(fp)) != EOF) {
			if (c == '\t' || c == '\n' || c == '\r') {
				str[i] = 0;
				if (str[0] != '@') {
					k = kh_get(str, h, str);
					if (k != kh_end(h))
						bns->anns[kh_val(h, k)].is_alt = 1;
				}
				while (c != '\n' && c != EOF) c = fgetc(fp);
				i = 0;
			} else str[i++] = c; // FIXME: potential segfault here
		}
		kh_destroy(str, h);
		fclose(fp);
	}
	return bns;
}

void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) err_fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

int bns_pos2rid(const bntseq_t *bns, int64_t pos_f)
{
	int left, mid, right;
	if (pos_f >= bns->l_pac) return -1;
	left = 0; mid = 0; right = bns->n_seqs;
	while (left < right) { // binary search
		mid = (left + right) >> 1;
		if (pos_f >= bns->anns[mid].offset) {
			if (mid == bns->n_seqs - 1) break;
			if (pos_f < bns->anns[mid+1].offset) break; // bracketed
			left = mid + 1;
		} else right = mid;
	}
	return mid;
}

int bns_intv2rid(const bntseq_t *bns, int64_t rb, int64_t re)
{
	int is_rev, rid_b, rid_e;
	if (rb < bns->l_pac && re > bns->l_pac) return -2;
	assert(rb <= re);
	rid_b = bns_pos2rid(bns, bns_depos(bns, rb, &is_rev));
	rid_e = rb < re? bns_pos2rid(bns, bns_depos(bns, re - 1, &is_rev)) : rid_b;
	return rid_b == rid_e? rid_b : -1;
}

int bns_cnt_ambi(const bntseq_t *bns, int64_t pos_f, int len, int *ref_id)
{
	int left, mid, right, nn;
	if (ref_id) *ref_id = bns_pos2rid(bns, pos_f);
	left = 0; right = bns->n_holes; nn = 0;
	while (left < right) {
		mid = (left + right) >> 1;
		if (pos_f >= bns->ambs[mid].offset + bns->ambs[mid].len) left = mid + 1;
		else if (pos_f + len <= bns->ambs[mid].offset) right = mid;
		else { // overlap
			if (pos_f >= bns->ambs[mid].offset) {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pos_f + len?
					bns->ambs[mid].offset + bns->ambs[mid].len - pos_f : len;
			} else {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pos_f + len?
					bns->ambs[mid].len : len - (bns->ambs[mid].offset - pos_f);
			}
			break;
		}
	}
	return nn;
}

uint8_t *bns_get_seq(int64_t l_pac, const uint8_t *pac, int64_t beg, int64_t end, int64_t *len)
{
	uint8_t *seq = 0;
	if (end < beg) end ^= beg, beg ^= end, end ^= beg; // if end is smaller, swap
	if (end > l_pac<<1) end = l_pac<<1;
	if (beg < 0) beg = 0;
	if (beg >= l_pac || end <= l_pac) {
		int64_t k, l = 0;
		*len = end - beg;
		seq = (uint8_t*)malloc(end - beg);
		if (beg >= l_pac) { // reverse strand
			int64_t beg_f = (l_pac<<1) - 1 - end;
			int64_t end_f = (l_pac<<1) - 1 - beg;
			for (k = end_f; k > beg_f; --k)
				seq[l++] = 3 - _get_pac(pac, k);
		} else { // forward strand
			for (k = beg; k < end; ++k)
				seq[l++] = _get_pac(pac, k);
		}
	} else *len = 0; // if bridging the forward-reverse boundary, return nothing
	return seq;
}

uint8_t *bns_fetch_seq(const bntseq_t *bns, const uint8_t *pac, int64_t *beg, int64_t mid, int64_t *end, int *rid)
{
	int64_t far_beg, far_end, len;
	int is_rev;
	uint8_t *seq;

	if (*end < *beg) *end ^= *beg, *beg ^= *end, *end ^= *beg; // if end is smaller, swap
	assert(*beg <= mid && mid < *end);
	*rid = bns_pos2rid(bns, bns_depos(bns, mid, &is_rev));
	far_beg = bns->anns[*rid].offset;
	far_end = far_beg + bns->anns[*rid].len;
	if (is_rev) { // flip to the reverse strand
		int64_t tmp = far_beg;
		far_beg = (bns->l_pac<<1) - far_end;
		far_end = (bns->l_pac<<1) - tmp;
	}
	*beg = *beg > far_beg? *beg : far_beg;
	*end = *end < far_end? *end : far_end;
	seq = bns_get_seq(bns->l_pac, pac, *beg, *end, &len);
	if (seq == 0 || *end - *beg != len) {
		fprintf(stderr, "[E::%s] begin=%ld, mid=%ld, end=%ld, len=%ld, seq=%p, rid=%d, far_beg=%ld, far_end=%ld\n",
				__func__, (long)*beg, (long)mid, (long)*end, (long)len, seq, *rid, (long)far_beg, (long)far_end);
	}
	assert(seq && *end - *beg == len); // assertion failure should never happen
	return seq;
}
