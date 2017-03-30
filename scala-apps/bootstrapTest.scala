if (args.size < 4 ) println("scala test #ref #alt #refObs #altObs")

val ran = scala.util.Random
val iterations = 10000

val dataPool = Array.fill[String](args(0).toInt)("R") ++ Array.fill[String](args(1).toInt)("A")
val reads = args(2).toInt + args(3).toInt
val tAR = args(3).toFloat / (args(2).toFloat + args(3).toFloat)
var bigOrEq, smallOrEq = 0
var tmp = 0
while (tmp < iterations){
//for (i <- 1 to iterations){
	var ref = 0
	var alt = 0
	var ctr = 0
	//for ( n <- 1 to reads){
	while (ctr < reads){
		if (dataPool(ran.nextInt(dataPool.size)) == "R") ref += 1 else alt += 1
		ctr += 1
	}
	if ((alt/(ref+alt).toDouble) < tAR) smallOrEq += 1 else if ((alt/(ref+alt).toDouble) >= tAR) bigOrEq += 1

	//results = (alt/(ref+alt).toDouble) :: results
	tmp += 1
}

println(smallOrEq)
println(bigOrEq)
println((if (smallOrEq == 0) 1/iterations.toDouble else smallOrEq/iterations.toDouble )+ " " + (if (bigOrEq == 0) 1/iterations.toDouble else bigOrEq/iterations.toDouble))
//println((if ((args(1).toInt/(args(0).toInt + args(1).toInt).toDouble) > 0.01) smallOrEq/iterations.toDouble else bigOrEq/iterations.toDouble))

/*
 + " " + 1/((args(0).toFloat + args(1).toFloat)/20)

val resultsSorted = results.toArray.sorted

println("Min = " + resultsSorted(0) + " Max = " + resultsSorted(iterations-1))

if (tAR < resultsSorted(0)){
	println("95% CI " + resultsSorted((iterations * 0.025).toInt) + "-" + resultsSorted((iterations*0.975).toInt) + " p < " + 1/iterations.toDouble)
} else{
	var count = 0
	while (resultsSorted(count) <= tAR && count < iterations-1) count += 1
	println("95% CI " + resultsSorted((iterations * 0.025).toInt) + "-" + resultsSorted((iterations*0.975).toInt) + " p = " + ((iterations-count)/iterations.toDouble))
}
*/