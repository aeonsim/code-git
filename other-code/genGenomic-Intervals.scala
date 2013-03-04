 import scala.io._
 import java.io._
 
 
 val data = Source.fromFile("chr-size.txt").getLines.toArray
 val out = new BufferedWriter(new FileWriter("Genomic_Intervals.txt"))
 
 val stepSize = 500000
 
 for (x <- data){
 val chrSize = x.split("\t")(1).toInt
 val chr = x.split("\t")(0)
 out.write("BASH_ARRAY0\n")
 var working = 0
 while (working < chrSize){
  if (working == 0){
  out.write(chr + ":" + (working + 1) + "-" + (working + stepSize) + "\n")
  } else {
 if ((working + stepSize) > chrSize){
 out.write(chr + ":" + (working + 1 - 5000) + "-" + chrSize + "\n")
 } else {
  out.write(chr + ":" + (working + 1 - 5000) + "-" + (working + stepSize) + "\n")
 }
 }
 working += stepSize
 }
 }
 out.close
 