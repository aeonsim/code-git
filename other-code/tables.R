setwd("/Users/chhar0/Documents/Phd-work/Code")
auto = read.table("autozygosity.prob")

setwd("/Users/chhar0/Documents/Phd-work/def")
data = read.table("file_of_interest.txt")
attach (data)
data = data[order(-V6)]
detach (data)

s1mb = subset(data, V6 > 1000000)
s5mb = subset(data, V6 > 5000000)
