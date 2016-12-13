import java.io._
if (args.size < 1 ) println("scala test input.txt ")
val ran = scala.util.Random

val input = new BufferedReader(new FileReader(new File(args(0))))

while (input.ready){

	val line = input.readLine.split("\t")
	System.err.println(line(1) + " " + line(2))

	val phase = line(6)

	val testCol = if (args.contains("proband")) 5 else if (phase == "M") 4 else 3


	if (line(testCol) != "-" && line(5) != "-"){

		val delim = if (line(testCol).contains(":")) ":" else ","

		val test = line(testCol).split(delim).map(_.toInt).sum 

		val tmpTest = line(testCol).split(delim).map(_.toDouble)

		val AR = tmpTest(1)/(tmpTest(0) + tmpTest(1))
		val testSize = tmpTest(0) + tmpTest(1)
		//System.err.println(testSize)

		var count = if (args.contains("proband")) 7 else 5
		var altAll, refAll = 0

			while (count < (line.size - 1)){
				if (phase == line(count + 1) && line(count) != "-"){
					val tmp = line(count).split(delim)
					if (tmp(1).toDouble/(tmp(1).toDouble + tmp(0).toDouble) >= 0.075){
						refAll += tmp(0).toInt
						altAll += tmp(1).toInt
					}
				}
				count += 2
			}


			var pool : List[String] = Nil

			for (i <- 1 to refAll) pool = "R" :: pool
			for (i <- 1 to altAll) pool = "A" :: pool

			val dataPool = if (pool.size == 0) Array("A","R") else pool.toArray
			var results : List[Double] = Nil
			//System.err.println("Datapool = " + dataPool.size)
			for (i <- 1 to 100000){
				var ref = 0
				var alt = 0
				for ( notUsed <- 1 to testSize.toInt){
					if (dataPool(ran.nextInt(dataPool.size)) == "R") ref += 1 else alt += 1
				}
				results = (alt/(ref+alt).toDouble) :: results
			}

			val resultsSorted = results.sorted

			if (AR < resultsSorted(0)){
				println(line(0) + "\t" + line(1) + "\t" + line(2) + "\t95%\t" + resultsSorted(2499) + "-" + resultsSorted(97499) + "\tp<\t1E-05\t" + refAll + "\t" + altAll + "\t" + AR + "\t" + tmpTest(0) + "\t" + tmpTest(1))
			} else{
				var count = 0
				while (count < 100000 && resultsSorted(count) <= AR) count += 1
				println(line(0) + "\t" + line(1) + "\t" + line(2) + "\t95%\t" + resultsSorted(2499) + "-" + resultsSorted(97499) + "\tp=\t" + (count/100000.0) + "\t" + refAll + "\t" + altAll + "\t" + AR + "\t" + tmpTest(0) + "\t" + tmpTest(1))
			}
	}

}