import java.io._
import scala.io._

object nucD {

def main(args:Array[String]) : Unit = {

for (num <- (1 to args(1).toInt)){

val data = new BufferedReader(new FileReader("autozygosity.prob"))
val dataOut = new BufferedWriter(new FileWriter("autozygous_blocks"+ num + ".txt"))

//val num = 23
var dataLine = data.readLine.trim.split(" +")(num).toDouble
var window = new Array[Double](100000)
var basePos = 0
var start = -1
val windowSize = 15

while (data.ready){
var count = 0
while ((count < 100000) && data.ready){
window(count) = dataLine
dataLine = data.readLine.trim.split(" +")(num).toDouble
count += 1
}
val current = window.iterator.sliding(windowSize).toArray
for (wind <- current){
if ((wind.sum / windowSize) >= 0.6){
if (start == -1) {start = current.indexOf(wind) + basePos}
} else {
if (start != -1){
dataOut.write(start + "\t" + (current.indexOf(wind)+basePos) + "\t" + (wind.sum / windowSize) + "\n" )
start = -1
}
}
}
basePos += count
}

dataOut.close
data.close

val dataIn = new BufferedReader(new FileReader("autozygous_blocks"+ num + ".txt"))
val idx = new BufferedReader(new FileReader("snp_indx.txt"))
val bedOut = new BufferedWriter(new FileWriter("out" + num + ".txt"))
val snpsIn = new BufferedReader(new FileReader("homozygosity_id.txt"))

import scala.collection.mutable.HashMap
var index = new HashMap[Int,Tuple2[String,Int]]
var snps = new HashMap[Int,Float]

while (snpsIn.ready){
val cline = snpsIn.readLine.trim.split(" +")
snps += cline(0).toInt -> cline(num).toFloat
}

while (idx.ready){
val cline = idx.readLine.split("\t")
index += cline(0).toInt -> (cline(1),cline(2).toInt)
}

while (dataIn.ready){
val cline = dataIn.readLine.split("\t")
if ((index.contains(cline(0).toInt)) && (index(cline(0).toInt)._1.toUpper != "CHRX")){
var het,hom = 0
var pos = cline(0).toInt
while (pos <= cline(1).toInt){
snps(pos) match {
case 0.0 => het += 1
case 1.0 => hom += 1
case _ => ;
}
pos += 1
}
bedOut.write(index(cline(0).toInt)._1 + "\t" + (index(cline(0).toInt)._2 - 1) + "\t" + (index(cline(1).toInt)._2 - 1) + "\t" + (het.toFloat/(het.toFloat + hom.toFloat)) + "\t" + (index(cline(1).toInt)._2 - index(cline(0).toInt)._2) + "\n")
}
}
bedOut.close
}

}
}