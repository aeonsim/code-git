//Find Mosaic dist
import java.io._
import scala.collection.mutable.HashMap

println("scala findMosaicDist.scala /path/to/gametes.tab minAR maxAR numberOfevents")

var in = new BufferedReader(new FileReader(new File(args(0))))

val min = args(1).toDouble
val max = args(2).toDouble
val number = args(3).toInt

var results = new HashMap[Int, Int]
var same = new HashMap[Int, Int]
var goodCand : List[String] = Nil

while (in.ready){
	val line = in.readLine.split("\t")
	var counter = 2
	var qual = 0
	var shared = 0

	while (counter <= line.size){
		if (line(counter).toDouble >= min && line(counter).toDouble <= max){
			qual += 1
			if (goodCand != Nil){
				if (goodCand.contains(line(counter - 1))) shared += 1
			}
		}
		counter += 2
	}

	if (qual == number && goodCand == Nil){
		counter = 2 
		while (counter <= line.size){
			goodCand = line(counter - 1) :: goodCand
			counter += 2
		}
	}

	if (results.contains(qual)) results(qual) += 1 else results += qual -> 1
	if (same.contains(shared)) same(shared) += 1 else same += shared -> 1


}

//println("Discovery")
//println("Dist of Mosaic variants")
//results.toArray.sorted.foreach(a => println(a._1 + "\t" + a._2))

//println("\nDist of Shared De novos")
//same.toArray.sorted.foreach(a => println(a._1 + "\t" + a._2))

//println("True values")

in = new BufferedReader(new FileReader(new File(args(0))))

results = new HashMap[Int, Int]
same = new HashMap[Int, Int]
var numLines = 0.0

while (in.ready){
	val line = in.readLine.split("\t")
	var counter = 2
	var qual = 0
	var shared = 0
	numLines += 1

	while (counter <= line.size){
		if (line(counter).toDouble >= min && line(counter).toDouble <= max){
			qual += 1
				if (goodCand.contains(line(counter - 1))) shared += 1
		}
		counter += 2
	}

	if (results.contains(qual)) results(qual) += 1 else results += qual -> 1
	if (same.contains(shared)) same(shared) += 1 else same += shared -> 1


}

println("Dist of Mosaic variants")
results.toArray.sorted.foreach(a => println("Dist\t" + a._1 + "\t" + a._2 + "\t" + a._2/numLines))

println("\nDist of Shared De novos")
same.toArray.sorted.foreach(a => println("Shared\t" + a._1 + "\t" + a._2 + "\t" + a._2/numLines))
println("end")
in.close