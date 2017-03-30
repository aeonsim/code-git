import java.io._
import scala.collection.mutable.HashSet

val inData = new BufferedReader(new FileReader(new File("freebayes.ab.tab")))
val keep = new BufferedReader(new FileReader(new File("ofInterest.txt")))

var keepSet = new HashSet[String]
var done = new HashSet[String]

while (keep.ready) keepSet += keep.readLine

var tmp = inData.readLine
println(tmp)

while (inData.ready){
	tmp = inData.readLine
	val tmpA = tmp.split("\t")
	if (keepSet.contains(tmpA(2) && ! done.contains(tmpA(2)))){
		done += tmpA(2)
		println(tmp)
	}
}


/* VCF filter */

import java.io._
import scala.collection.mutable.HashSet
import net.sf.samtools.util._

val inData = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val keep = new BufferedReader(new FileReader(new File(args(1))))

var keepSet = new HashSet[String]
//var done = new HashSet[String]

while (keep.ready) keepSet += keep.readLine

var tmp = inData.readLine
while(tmp(0) == '#'){
	println(tmp)
	tmp = inData.readLine
}


while (inData.ready){
	tmp = inData.readLine
	val tmpA = tmp.split("\t")
	if (keepSet.contains(tmpA(0) + ":" + tmpA(1))){
		//done += tmpA(2)
		println(tmp)
	}
}