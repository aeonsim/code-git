object gcD{

	def tstv (ts: Array[String], tv: Array[String], events: List[String]): Unit = {
		var counter = 0
		var tsNum, tvNum = 0.0
		while (counter < events.size){
			if (ts.contains(events(counter))) tsNum += 1.0 else if (tv.contains(events(counter))) tvNum += 1.0
			counter += 1
		}
		println(" TS/TV is: " + tsNum/tvNum)

	}


	def main (args: Array[String]) : Unit = {
		import java.io._

		val tsA = Array("AG","GA","CT","TC")
		val tvA = Array("CG","GC","AT","TA","AC","CA","GT","TG")
		val generations = 100000
		var snps : List[String] = Nil
		var denovos : List[String] = Nil

		val in_denovos = new BufferedReader(new FileReader(new File(args(0))))
		val in_snps = new BufferedReader(new FileReader(new File(args(1))))

		val gcBias = 67
		val gcEvents = 158*5
		val denovoEvents = 50*158/2600

		while (in_snps.ready){
			var snpLine = in_snps.readLine.split("\t")
			while (snpLine(0)(0) == '#') in_snps.readLine.split("\t")
			if ( snpLine(3).size == 1 && snpLine(4).size == 1 ){
				snps = (snpLine(3) + snpLine(4)) :: snps
			}
		}

		in_snps.close

		while (in_denovos.ready){
			denovos = in_denovos.readLine :: denovos
		}

		in_denovos.close

		var cgen = 0

		while (cgen < generations ){
			print("Generation: " + cgen)
			tstv(tsA,tvA,snps)
			for ( i <- 1 to denovoEvents){
				snps = denovos(scala.util.Random.nextInt(denovos.size)) :: snps
			}

			

		}



	}
}