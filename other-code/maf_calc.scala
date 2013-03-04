import scala.io.Source._
import java.io._

val in = fromFile("genos.txt").getLines.toArray
val out = new BufferedWriter(new FileWriter("mafs.txt"))
val GTpos = 0

var x = 0
 while (x < in.size){
 var AA, BB = 0
 val cline = in(x).split("\t")
 for (anml <- cline){
 var anlist = anml.split(":")
 anlist(GTpos) match {
 case "0/0" => AA += 2
 case "0/1" => AA += 1; BB += 1
 case "1/0" => AA += 1; BB += 1 
 case "1/1" => BB += 2
 case "./." =>; 
 case _ =>;
 }
 }
 val MAF = if ((AA + BB) == 0) 0 else AA.toFloat / (AA + BB).toFloat
out.write((x + 1) + "  " + MAF + "\n")
 x += 1
 }
 out.close
 
 
 
 /*
var x = 0
while (x < in.size){

var AA, BB = 0
for (cur <- in(x).split("\t")){
cur match {
case "1" => AA += 2
case "2" => AA += 1; BB += 1
case "3" => BB += 2
case _ =>;
}
}
val MAF = if ((AA + BB) == 0) 0 else AA.toFloat / (AA + BB).toFloat
out.write((x + 1) + "  " + MAF + "\n")
//out.write(MAF + "\n")
x += 1
}
out.close

*/