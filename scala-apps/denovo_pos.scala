val chr1 = 158337067

import java.util.Random

val ran = new Random
var avgs: List[Int] = Nil

println("Dif\tGroup")

for (y <- 0 to 10000){
var tmpList: List[Int] = Nil
for (i <- 0 to 4){
tmpList = ran.nextInt(158337067) :: tmpList
}
val tmp = tmpList.sorted
for (i <- 1 to 4){
avgs = (tmp(i) - tmp(i -1)) :: avgs
}
}

avgs.foreach(s => println(s + "\tsim"))
