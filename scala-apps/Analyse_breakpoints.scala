import java.io._
import scala.collection.mutable.HashMap

val in = new BufferedReader(new FileReader(new File("/Users/chhar0/tmp/input.reads.txt")))
var data = new HashMap[String,List[String]]

val outLeft = new BufferedWriter(new FileWriter(new File("ERV_5p.tab")))
val outRight = new BufferedWriter(new FileWriter(new File("ERV_3p.tab")))

while (in.ready){
	val tmp = in.readLine.split("\t")
	val flags = tmp(0)
	val cigar = tmp(1)
	val read = tmp(2)
	val id = tmp(3).split(":").apply(2)
	if (cigar.contains("S") && cigar.contains("M") && cigar.map(s => if (s == 'M') 1 else 0).sum == 1){

		if(cigar.indexOf("S") < cigar.indexOf("M")){
			//case Softclip 5'
			val bases = cigar.substring(0,cigar.indexOf("S")).toInt
			if (bases >= 10){
				outLeft.write(read.substring(bases -10, bases) + "__" + read.substring(bases,bases + 10) + "\t" + id + "\t" + flags + "\n")
			} else {
				for (i <- 0 to (scala.math.abs(bases - 10) -1)) outLeft.write("N")
				outLeft.write(read.substring(0, bases) + "__" + read.substring(bases,bases + 10) + "\t" + id + "\t" + flags + "\n")
			}

		} else {
				//case softclip 3'
			val bases = cigar.substring(0,cigar.indexOf("M")).toInt
			val split = cigar.substring(cigar.indexOf("M")+1,cigar.indexOf("S")).toInt
			if (split >= 10){
				outRight.write(read.substring(bases -10, bases) + "__" + read.substring(bases,bases + 10) + "\t" + id + "\t" + flags + "\n")
			} else {
				outRight.write(read.substring(bases - 10, bases) + "__" + read.substring(bases,bases + split))
				for (i <- 0 to (scala.math.abs(split - 10) -1)) outRight.write("N")
				outRight.write("\t" + id + "\t" + flags + "\n")
			}

		}

	}

}