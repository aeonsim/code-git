import scala.io.Source._
import java.io._

object tobed {

def main( args: Array[String]){

val anidx = fromFile("anml_indx.txt").getLines.toArray

for (anml <- (0 to anidx.size -1).par){ 

val name = anidx(anml - 1).split("\t").apply(1)
val aut = new BufferedReader(new FileReader("autozygosity.prob"))
val idx = new BufferedReader(new FileReader("snp_indx.txt"))
val out = new BufferedWriter(new FileWriter("animal_" + name + ".bedgraph"))

out.write("track type=bedGraph name=Autozyg_" + name + "description=\"HAutozygous Prob\" autoscale=off viewLimits=0:1\n")

while(aut.ready){
val auto = aut.readLine.trim.split(" +")
val snps = idx.readLine.split("\t")
//println(auto(anml)+ " " + auto(1)+ " " + auto(0) + " " + auto(2))
out.write(snps(1) + "\t" + (snps(2).toInt - 1) + "\t" + snps(2) + "\t" + auto(anml) + "\n")
}
aut.close
idx.close
out.close
}

}
}