-- Proteins Table
create table proteins (`protein_id` String , `size` UInt64, `helices` UInt64, `sheets` UInt64,`loops` UInt64,`bends` UInt64,`sequence` String,`ss` String, `pattern` String,`enhanced` String,`processdate` Date DEFAULT today()) Engine=MergeTree(processdate,(protein_id,helices,sheets,loops,bends,sequence,ss,pattern,enhanced),8192);
-- Helix
create table structures (`protein_id` String , `ss_id` String,`type` String,`pos` UInt64, `size` UInt64,`sequence` String, `residues` Array(UInt8),`sasa` Array(Float32),`pattern` Array(UInt8),`processdate` Date DEFAULT today()) Engine=MergeTree(processdate,(protein_id,ss_id,type,pos,size,sequence,residues,sasa,pattern),8192);
-- -- sheets
-- create table sheets (`protein_id` String , `ss_id` String, `size` UInt64, `residues` Array(UInt64),`pattern` Array(UInt64),`processdate` Date DEFAULT today()) Engine=MergeTree(processdate,(protein_id,ss_id,size,residues,pattern),8192);
-- -- loops
-- create table sheets (`protein_id` String , `ss_id` String, `size` UInt64, `residues` Array(UInt64),`pattern` Array(UInt64),`processdate` Date DEFAULT today()) Engine=MergeTree(processdate,(protein_id,ss_id,size,residues,pattern),8192);
-- --bends
-- create table sheets (`protein_id` String , `ss_id` String, `size` UInt64, `residues` Array(UInt64),`pattern` Array(UInt64),`processdate` Date DEFAULT today()) Engine=MergeTree(processdate,(protein_id,ss_id,size,residues,pattern),8192);

