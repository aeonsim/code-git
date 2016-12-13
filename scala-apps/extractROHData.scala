import org.apache.commons.io.FileUtils._
import scala.collection.immutable.HashMap
import java.io._

val fwd = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("txt"),true).iterator

var first = true
var pos: List[String] = Nil
var data : List[Array[String]] = Nil
var fpos = new HashMap[String,Int]

while (fwd.hasNext){

	val F = fwd.next

	val fID = F.toString.split("/").last.split("_").last.split("\\.").apply(0)
	var curdata : List[String] = Nil
	val in = new BufferedReader(new FileReader(F))
	var ctmp = in.readLine
	while (ctmp(0) == '#'){
			ctmp = in.readLine
			//println(ctmp)
		} 

	if (first){
		curdata = fID :: curdata
		while (in.ready){
			val tmp = in.readLine.split("\t")
			//println(tmp)
			pos = tmp(1) :: pos
			curdata = tmp(2) :: curdata
		}
		System.err.println("first done")
		for (i <- pos) fpos += i -> 0
		first = false
	} else {
		curdata = fID :: curdata
		while (in.ready){
			val tmp = in.readLine.split("\t")
			if (fpos.contains(tmp(1))) curdata = tmp(2) :: curdata
		}

	}

	data = curdata.reverse.toArray :: data
	System.err.println(fID + " loaded " + data.size)
}

var count = 0
var orderedPos = pos.reverse.toArray

System.err.println(data.size)
System.err.println(data(0).size)

while (count < orderedPos.size){
if (count == 0) print("...\t") else print(orderedPos(count) + "\t")
for (i <- data){
//System.err.println(i.contains(orderedPos(count)) + " pos is " + orderedPos(count))
print(i(count) + "\t")
}
print("\n")
count += 1
}