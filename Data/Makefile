DAT := $(shell ls -1 | awk '/^CRR/')
HDR := $(foreach x, $(DAT), processed/$(x)/counts_filtered)

EXT := mtx.gz rows.gz cols.gz
FULL := $(foreach hdr, $(HDR), $(foreach x, $(EXT), $(hdr).$(x)))

all: $(FULL)

processed/%.mtx.gz: swap/%/spliced.mtx.gz swap/%/spliced.mtx.gz.index swap/%/spliced.rows.gz swap/%/spliced.cols.gz swap/%/unspliced.mtx.gz swap/%/unspliced.mtx.gz.index swap/%/unspliced.rows.gz swap/%/unspliced.cols.gz
	mkdir -p $(dir $@)
	Rscript --vanilla -e "mmutilR::vcat.sparse('swap/$*/spliced', 'swap/$*/unspliced', 'processed/$*')"

processed/%.rows.gz: processed/%.mtx.gz
	mkdir -p $(dir $@)
	touch $@

processed/%.cols.gz: processed/%.mtx.gz
	mkdir -p $(dir $@)
	touch $@

## Since Kallisto-BUStools generate cell x gene mtx files, we would
## need to swap row and column position in the KB-generated mtx
swap/%.mtx.gz: %.mtx
	mkdir -p $(dir $@)
	cat $^ | awk '/^%/ { print } !/^%/{ print $$2 FS $$1 FS $$3 }' | bgzip -c > $@

swap/%.mtx.gz.index: swap/%.mtx.gz
	echo $@

swap/%.rows.gz: %.genes.txt
	mkdir -p $(dir $@)
	cat $< | awk -vTAG=$$(basename $@ .rows.gz) '{ print $$0 "_" TAG }' | gzip -c > $@

swap/%.cols.gz: %.barcodes.txt
	mkdir -p $(dir $@)
	cat $< | gzip -c > $@

