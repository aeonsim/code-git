import java.io._

System.err.println("scala app data.txt name")

val in = new BufferedReader(new FileReader(new File(args(0))))

var data : List[Double] = Nil
var names : List[String] = Nil

while (in.ready){
	val tmp = in.readLine.split("\t")
	data = tmp(1).toDouble :: data
	names = tmp(0) :: names
}

var proAvg, offAvg = 0.0
var proNum, offNum = 0

for (i <- 0 to (data.size -1)){
	if (names(i) == args(1)){
		proAvg += data(i)
		proNum += 1
	} else {
		offAvg += data(i)
		offNum += 1
	}
}

val obsDif = (proAvg / proNum) - ( offAvg / offNum )

println(s"${args(1)} mean AR is: ${proAvg / proNum}\nOffspring mean AR is: ${offAvg / offNum}\nObserved difference ${obsDif}")

var count = 0

var proSimAvg, offSimAvg, avgDif : List[Double] = Nil

while (count < 1000000){
	val ldata = scala.util.Random.shuffle(data)
	val lnames = scala.util.Random.shuffle(names)

	var lproAvg, loffAvg = 0.0
	var lproNum, loffNum = 0

	for (i <- 0 to (ldata.size -1)){
		if (lnames(i) == args(1)){
			lproAvg += ldata(i)
			lproNum += 1
		} else {
			loffAvg += ldata(i)
			loffNum += 1
		}
	}

	avgDif = (lproAvg / lproNum) - ( loffAvg / loffNum ) :: avgDif
	count += 1
}

val sortedDif = avgDif.sorted

if (obsDif < sortedDif(0)){
	println("95%\t" + sortedDif(24999) + "-" + sortedDif(974999) + "\tp<\t1E-06\t" + " Min: " + sortedDif(0) + "\tMax: " + sortedDif(999999))
} else{
	var count2 = 0
	while (count2 < 1000000 && sortedDif(count2) <= obsDif) count2 += 1
	println("95%\t" + sortedDif(24999) + "-" + sortedDif(974999) + "\tp=\t" + (count2/1000000.0) + " Min: " + sortedDif(0) + "\tMax: " + sortedDif(999999))
}