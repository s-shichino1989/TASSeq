import hts
import os
import strutils
import tables
var ibam:Bam

# lookup from cb -> cluster
var clusterTbl = initTable[string,string]()
# lookup from cluster -> bam
var tbl = initTable[string, Bam]()

for x in paramStr(1).lines:
  var toks = x.strip().split(",")
  clusterTbl[toks[0]] = toks[1]

if not open(ibam, paramStr(2)):
   quit "couldn't open bam"

var filename_prefix = paramStr(3)

for aln in ibam:
  var cb = tag[string](aln, "CB").get
  if cb.isEmptyOrWhitespace: continue
  if cb notin clusterTbl: continue
  var cluster = clusterTbl[cb]
  if cluster notin tbl:
    var obam: Bam
    if not open(obam, cluster & ".bam", mode="w", threads=4):
      quit "couldn't open bam for writing"
    obam.write_header(ibam.hdr)
    tbl[cluster] = obam
  tbl[cluster].write(aln)

for k, bam in tbl:
  bam.close()
ibam.close()
