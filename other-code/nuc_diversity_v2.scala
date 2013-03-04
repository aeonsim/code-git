import java.io._
import scala.io._
import scala.collection.mutable.HashMap

object nucD {

def main(args:Array[String]) : Unit = {

for (num <- (1 to args(0).toInt)){

val data = new BufferedReader(new FileReader("autozygosity.prob"))
//val dataOut = new BufferedWriter(new FileWriter("autozygous_blocks"+ num + ".txt"))
var autozyg : List[Tuple2[Int,Int]] = Nil

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
autozyg = (start, current.indexOf(wind)+basePos) :: autozyg
//dataOut.write(start + "\t" + (current.indexOf(wind)+basePos) + "\t" + (wind.sum / windowSize) + "\n" )
start = -1
}
}
}
basePos += count
}

//dataOut.close
data.close

//val dataIn = new BufferedReader(new FileReader("autozygous_blocks"+ num + ".txt"))
val idx = new BufferedReader(new FileReader("snp_indx.txt"))
val bedOut = new BufferedWriter(new FileWriter("out" + num + ".txt"))
val snpsIn = new BufferedReader(new FileReader("homozygosity_id.txt"))

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

for (auto <- autozyg.reverse){
if ((index.contains(auto._1)) && (index(auto._1)._1.toUpperCase != "CHRX")){
var het,hom = 0
var pos = auto._1
while ((pos <= auto._2)&&(snps.contains(pos))){
snps(pos) match {
case 0.0 => het += 1
case 1.0 => hom += 1
case _ => ;
}
pos += 1
}
if ((index.contains(auto._1))&&(index.contains(auto._2))&&(index(auto._1)._1 == index(auto._2)._1) {
bedOut.write(index(auto._1)._1 + "\t" + (index(auto._1)._2 - 1) + "\t" + (index(auto._2)._2 - 1) + "\t" + (het.toFloat/(het.toFloat + hom.toFloat)) + "\t" + het + "," + hom + "\t" + (index(auto._2)._2 - index(auto._1)._2) + "\n")
} else {
bedOut.write(index(auto._1)._1 + "\t" + (index(auto._1)._2 - 1) + "\t" + (index(auto._2 - 1)._2 - 1) + "\t" + (het.toFloat/(het.toFloat + hom.toFloat)) + "\t" + het + "," + hom + "\t" + (index(auto._2)._2 - index(auto._2 - 1)._2) + "\n")
}
}
}
bedOut.close
}
}
}