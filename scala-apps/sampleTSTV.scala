import java.io._

val in = new BufferedReader(new FileReader(new File(args(0))))
val sample = args(1).toInt
val testTsTv = args(2).toDouble
val trials = 100000

var tmpData : List[String] = Nil

while (in.ready){
	val tmp = in.readLine.split("\t")
	tmpData = tmp(0) + ">" +  tmp(1) :: tmpData
}

in.close
val data = tmpData.toArray

var count = 0

val ts = Array("C>T","G>A","T>C","A>G")
val tv = Array("C>A","A>C","C>G","G>C","A>T","T>A","T>G","G>T")

var results : List[Double] = Nil

while (count < trials){
	count += 1
	var tsNum, tvNum, mutCount = 0
	while (mutCount < sample){
		mutCount += 1
		val mut = data(scala.util.Random.nextInt(data.size))
		if (ts.contains(mut)) tsNum += 1 else if (tv.contains(mut)) tvNum += 1
	}
	results = tsNum/tvNum.toDouble :: results
	System.err.println(tsNum/tvNum.toDouble)
}

val lessThan = results.filter(s => s <= testTsTv)
println(s"Number of Trials less than or Equal to ${testTsTv} is:\t${lessThan.size}\t${lessThan.size/trials.toDouble}")

var window = 0.0

while (window <= 4.5){
	window += 0.5
	val lt = results.filter(s => s <= window)
	println(s"${window}\t${lt.size}")
}