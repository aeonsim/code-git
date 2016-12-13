object chrComp {

	def main(args:Array[String]) : Unit = {
		System.err.println("scala app.scala numMutations mutationList.txt")
		val Chroms = Array(158337067,137060424,121430405,120829699,121191424,119458736,112638659,113384836,105708250,104305016,107310763,91163125,84240350,84648390,85296676,81724687,75158596,66004023,64057457,72042655,71599096,61435874,52530062,62714930,42904170,51681464,45407902,46312546,51505224,148823899)
		var corChrs : List[Int] = Nil

		import java.io._

		val input = new BufferedReader(new FileReader(new File(args(0))))
		var obsDataDis =  Array.fill[List[Int]](30)(Nil)
		var obsData = Array.fill[List[Int]](30)(Nil)
		var numMuts = 0
		while (input.ready){
			val line = input.readLine.split("\t")
			val chr = if(line(0).toUpperCase.contains("X")) 29 else if (line(0).toUpperCase.contains("CHR")) line(0).substring(3,line(0).size).toInt -1 else -1
			obsData(chr) = line(1).toInt :: obsData(chr)
			numMuts += 1
		}
		input.close
		println(numMuts)

		val numChrs = obsData.filter(s => s.size >= 1).size
		var mutChrs : List[Int] = Nil

		var i = 0
		while (i < Chroms.size){
			val chances = Chroms(i)/10000000
			for (p <- 1 to chances) corChrs = i :: corChrs
			i += 1
		}

		var dataDis = Array.fill[List[Int]](30)(Nil)
		var dataNum = Array.fill[List[Int]](30)(Nil)
		val rand = scala.util.Random
		var x = 0

		while (x < 1000000){
			var cRun = Array.fill[List[Int]](30)(Nil)
			var y = 0
			while (y <= numMuts){
				var cChr = corChrs(rand.nextInt(corChrs.size))
				cRun(cChr) = rand.nextInt(Chroms(cChr)) :: cRun(cChr)
				y += 1
			}
			//println("data gen")
			y = 0
			while (y < Chroms.size){
				if (cRun(y).size >= 2){
						var k = 0
						val ordered = cRun(y).sorted
						while (k < (ordered.size - 1)){
							dataDis(y) = (ordered(k+1) - ordered(k)) :: dataDis(y)
							k += 1
						}
						dataNum(y) = cRun(y).size :: dataNum(y)
					} else {
						dataNum(y) = cRun(y).size :: dataNum(y)
					}
				y += 1
			}

			mutChrs = cRun.filter(s => s.size >= 1).size :: mutChrs

			x += 1
		}
		println("sims")

		var c = 0
		while (c < Chroms.size){
			if (obsData(c).size >= 2){
				var k = 0
				val ordered = obsData(c).sorted
				while (k < (ordered.size - 1)){
					obsDataDis(c) = (ordered(k+1) - ordered(k)) :: dataDis(c)
					k += 1
				}
			}
			c += 1
		}
		println("data processed")

		c = 0
/*
		while (c < Chroms.size){
			if (obsDataDis(c).size >= 1){
				val avgDis = obsDataDis(c).sum / obsDataDis(c).size.toDouble
				val moreExtreme = dataDis(c).filter(s => s <= avgDis)
				println(s"chr${c+1}\t${moreExtreme.size/dataDis.size.toDouble}\t${moreExtreme.size}/${dataDis(c).size}")
			} else {
				val moreExtreme = dataDis(c).filter(s => s <= 0)
				println(s"chr${c+1}\t${0.0}\t${moreExtreme.size}/${dataDis(c).size}")
			}

			c += 1
		}
*/

		c = 0
		while (c < Chroms.size){
			val obsMuts = obsData(c).size
			val moreExtreme = dataNum(c).filter(s => s <= obsMuts)
			val moreExtremeGreater = dataNum(c).filter(s => s >= obsMuts)
			println(s"chr${c+1}\t${moreExtreme.size/dataNum(c).size.toDouble}\t${moreExtreme.size}/${dataNum(c).size}\t${moreExtremeGreater.size/dataNum(c).size.toDouble}\t${moreExtremeGreater.size}/${dataNum(c).size}")
			c += 1
		}

		val moreEx = mutChrs.filter(s => s <= numChrs) // list of number of simulated chrs with mut
		val moreExGreater = mutChrs.filter(s => s >= numChrs)

		println(s"Number of Obs Chrs with Mutations\t${numChrs}\t${moreEx.size}/1000000\t ${moreEx.size/1000000.0}\t${moreExGreater.size}/1000000\t${moreExGreater.size/1000000.0}")
		//numChrs // number of observed chrs with mut

	}


}
