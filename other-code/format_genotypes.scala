import scala.io.Source._
import java.io._

val in = fromFile("genos.txt").getLines.toArray
val out = new BufferedWriter(new FileWriter("genos_formated.txt"))

val GTpos = 0

var x = 0
 while (x < in.size){
 val cline = in(x).split("\t")
 for (anml <- cline){
 
 var anlist = anml.split(":")
 anlist(GTpos) match {
 case "0/0" => out.write("1\t")
 case "0/1" => out.write("2\t")
 case "1/1" => out.write("3\t")
 case _ => out.write("-1\t")
 }
 }
 out.write("\n")
 x += 1
 }
 out.close