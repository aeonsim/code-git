import java.io._
import scala.collection.mutable.HashMap

val in = new BufferedReader(new FileReader(new File(args(0))))

val R1 = new BufferedWriter(new FileWriter(new File("R1.fa")))
val R2 = new BufferedWriter(new FileWriter(new File("R2.fa")))

var data = new HashMap[String,HashMap[Int,Int]]

for (i <- List("chr1","chrX")) {
data += i -> HashMap(0 -> 0)
for (j <- 1 to 50) data(i) += (40 * j) -> 0
}

var cRead = 0

while (in.ready){
	val tmp = in.readLine.split("\t")
	val read = tmp(9)
	val length = (read.size /40).toInt * 40
	if (List("chr1","chrX").contains(tmp(2)) && ! tmp(5).contains("L")){
		if (length > 2000) data(tmp(2)).apply(2000) += 1 else data(tmp(2)).apply(length) += 1
		val reads = read.splitAt(read.length/2)

		R1.write(s">Read_${cRead}_R/1\n${reads._1}\n")
		R2.write(s">Read_${cRead}_R/2\n${reads._2}\n")	
	}
	
	cRead += 1
}

R1.close
R2.close
in.close
println("Size\tChr1\tChrX")
for ( i <- 0 to 50){
print(s"${i*40}\t${data("chr1")(i*40)}\t${data("chrX")(i*40)}\n")

}
