import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import org.apache.commons.io.FileUtils._
import java.io._


class PSEUDO(chr: String, startP: Int, endP: Int){
	var chrom = chr
	var start = startP
	var end = endP
	var fwdIP = 0
	var revIP = 0
	var fwdSR = 0
	var revSR = 0
	var mateSR = 0
	var difs : List[Int] = Nil
	var rPmP, rPmN, rNmP, rNmN = 0
	var events = new HashMap[String,Tuple2[Int,Int]] // Event, fwd, rev reads
	var geno = new HashMap[String,Tuple2[Int,Int]]  //indvidual, Depth, Bridging reads
	var windowDepth: List[Int] = Nil
	var other : List[Int] = Nil
	var etype : List[String] = Nil
	var best5p, best3p : List[Int] = Nil
	var size : List[Int] = Nil
	var readsRatio : List[Double] = Nil
	var carriers : List[String] = Nil
	var allBreaks = new HashMap[Int,Int]
	var genes = new HashSet[String]

	def mergeData(old: HashMap[String,Tuple2[Int,Int]], newD: HashMap[String,Tuple2[Int,Int]]): HashMap[String,Tuple2[Int,Int]] = {
		var report = old
		for (d <- newD){
			if (old.contains(d._1)){
				val tOld = old(d._1)
				val tNew = d._2

				report(d._1) = (tOld._1 + tNew._1,tOld._2 + tNew._2)
			} else {
				report += d._1 -> d._2
			}
		}
		report
	}

	def checkSame(chr: String, pos: List[Int]) : Boolean = {
		if (chr == chrom && ( pos.contains(start) || pos.contains(end))) true else false
	}

	def addRMInfo(vals: Array[Int]): Unit = {
		rPmP += vals(0)
		rPmN += vals(1)
		rNmP += vals(2)
		rNmN += vals(3)
	}

	def addEventInfo(vals : Array[Array[String]]): Unit = {
		for (x <- vals){
			if(events.contains(x(0))){
					events(x(0)) = (events(x(0))._1 + x(1).toInt, events(x(0))._2 + x(2).toInt)
				} else {
					events += x(0) -> (x(1).toInt,x(2).toInt)
				}
		}
	}

	def addGeno(id: String, data: String): Unit = {
		val tmp = data.split(":").map(_.toInt)
		geno += id -> (tmp(0),tmp(1))
	}

	def addBreaks(vals: String): Unit = {
		val tmp = vals.split(":").map(_.toInt)
		best5p = tmp(0) :: best5p
		best3p = tmp(1) :: best3p
	}

	def addSize(vals: String): Unit = {
		val tmp = if (vals != "NA") vals.split("-").map(_.toInt) else Array(-1,-1)
		if (tmp.size == 2) size = tmp(1) - tmp(0) :: size
	}

	def addMSR(num: String): Unit = {
		mateSR = mateSR + num.toInt
	}

}



object survey2 {

/*
	def convert2Pos(input: String): Array[Array[Int]] = {
		if (input != "Map()") input.slice(4,input.size -1).split(",").map(y => y.split("->").map(_.trim.toInt)) else Array()
	}
	
	def addBreaksInfo (breaks: HashMap[Int, Int], newBreaks: String ): HashMap[Int,Int] ={

		val novoBreaks = convert2Pos(newBreaks)// newBreaks.substring(4,newBreaks.size -1).split(",").map(s => s.trim.toInt)
		for (split <- novoBreaks){
			if (breaks.contains(split(0))){
				breaks(split(0)) = breaks(split(0)) + split(1)
			} else {
				breaks += split(0) -> split(1)		
			}
		}

		breaks
	}

	*/

	def convert2Pos(input: String): Array[Int] = {
		if (input != "Set()") input.slice(4,input.size -1).split(",").flatMap(y => y.split("->").map(_.trim.toInt)) else Array()
	}
	
	def addBreaksInfo (breaks: HashMap[Int, Int], newBreaks: String ): HashMap[Int,Int] = {

		val novoBreaks = convert2Pos(newBreaks)// newBreaks.substring(4,newBreaks.size -1).split(",").map(s => s.trim.toInt)
		for (split <- novoBreaks) {
			if (breaks.contains(split)) {
				breaks(split) = breaks(split) + 1
			} else {
				breaks += split -> 1		
			}
		}

		breaks
	}



	val expReads = 4035
	val expIPP = 88
	val expSR = 24

	def bootStrap(reads: Int, evid: Int, geno: Int, avg: Int): Double = {
		val ran = scala.util.Random
		val iterations = 10000
	    val dataPool = Array.fill[String](geno)("R") ++ Array.fill[String](avg)("A")
	    val tAR = evid.toDouble / reads.toDouble
	    var bigOrEq, smallOrEq = 0
	    var tmp = 0
	    if (dataPool.size >= 1) {
	      while (tmp < iterations) {
	        var ref = 0
	        var alt = 0
	        var ctr = 0
	        while (ctr < reads) {
	          if (dataPool(ran.nextInt(dataPool.size)) == "R") ref += 1 else alt += 1
	          ctr += 1
	        }
	        if (alt >= evid) bigOrEq += 1
	        tmp += 1
	      }

	      if (bigOrEq == 0) 1 / iterations.toDouble else bigOrEq / iterations.toDouble
	    } else {
	      1
	    }


	}



	def main(args: Array[String]) : Unit = {

		val dataSet = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("event.tab"),true).iterator
		val depths = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/currentDepths.txt")))
		val ped_in = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/CRV-pedigree-800K.ped")))
		val gap_in = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/gaps.bed")))

		val annoIn = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/Bos_taurus.UMD3.1.82.gtf")))
		val annoIn2 = new BufferedReader(new FileReader(new File("/home/u/u218757/kos_home/refs/RefSeq_Genes.gtf")))

		//Chr -> Tuple5[start, end, class, name/id, strand]

		var anno = new HashMap[String,Array[Tuple5[Int,Int,String,String,String]]]
		var annoList : List[Tuple5[Int,Int,String,String,String]] = Nil
		var prevChrom = ""
		
		while (annoIn.ready){
			val tmp = annoIn.readLine.split("\t")
			if (tmp(0).apply(0) != '#'){
				val chr = "chr" + tmp(0).trim
				if (tmp(2) == "transcript") {
					val info = tmp(8).split(" ")
					val name = if (info.contains("gene_name")) info(info.indexOf("gene_name") + 1 ) else info(info.indexOf("gene_id") + 1 )
					if (prevChrom == chr){
						annoList = (tmp(3).toInt,tmp(4).toInt, tmp(2), name, tmp(6)) :: annoList
					} else {
						anno += prevChrom -> annoList.reverse.toArray
						annoList = (tmp(3).toInt,tmp(4).toInt, tmp(2), name, tmp(6)) :: Nil
						prevChrom = chr
					}
				}
			}
		}//end while repeats
		
		anno += prevChrom -> annoList.reverse.toArray
		annoList = Nil
		prevChrom = ""
		annoIn.close
		
		
		
		var anno2 = new HashMap[String,Array[Tuple5[Int,Int,String,String,String]]]
		while (annoIn2.ready){
			val tmp = annoIn2.readLine.split("\t")
			if (tmp(0).apply(0) != '#'){
				val chr = tmp(0)
				if (tmp(2) == "gene") {
					val info = tmp(8).replaceAll("=",";").split(";")
					val name = if (info.contains("symbol_ncbi")) info(info.indexOf("symbol_ncbi") + 1) else if (info.contains("Dbxref")) info(info.indexOf("Dbxref") + 1) else info(info.indexOf("ID") + 1)
					if (prevChrom == chr){
						annoList = (tmp(3).toInt,tmp(4).toInt, tmp(2), name, tmp(6)) :: annoList
					} else {
						anno2 += prevChrom -> annoList.reverse.toArray
						annoList = (tmp(3).toInt,tmp(4).toInt, tmp(2), name, tmp(6)) :: Nil
						prevChrom = chr
					}
				}
			}
		}//end while repeats
		
		anno2 += prevChrom -> annoList.reverse.toArray
		annoList = Nil
		prevChrom = ""
		annoIn.close


//boolean, strand, Name?
		def isEle (pos: Int, curSearch: Array[Tuple5[Int,Int,String,String,String]]) : Tuple3[Boolean,String,String] = {
			val half = scala.math.floor(curSearch.length/2).toInt
			if (curSearch.length == 0){
					(false,"","")
			} else {
				if (curSearch(half)._1 <= pos && curSearch(half)._2 >= pos) {
					(true,curSearch(half)._5,s"${curSearch(half)._3}-${curSearch(half)._4}")
					} else {
						if (curSearch.length <= 1) {
							(false,"","")
						} else {
						if (curSearch(half)._1 > pos) isEle(pos,curSearch.slice(0,half)) else isEle(pos,curSearch.slice(half + 1,curSearch.length))
					}
				}
			}
		} //end isEle

		def isGene(pos: Int, chr: String): String = {
			val ebi = isEle(pos,anno(chr))
			val refSeq = isEle(pos,anno2(chr))
			if (ebi._1 == true){
					if (refSeq._1 == true){
						s"${ebi._2}|${ebi._3}|${refSeq._2}|${refSeq._3}"
					} else {
						s"${ebi._2}|${ebi._3}"
					}
				} else {
					if(refSeq._1 == true) s"${refSeq._2}|${refSeq._3}" else ""
				}
		}

		var ped = new HashMap[String,Array[String]]
		while (ped_in.ready){
			val tmp = ped_in.readLine.split("\t")
			ped += tmp(1) -> tmp
		}
		ped_in.close




		val gapsChrom = new HashMap[String,Array[Tuple2[Int,Int]]]
		var gapList: List[Tuple2[Int,Int]] = Nil
		//var prevChrom = ""

		while(gap_in.ready){
			val tmp = gap_in.readLine.split("\t")
			if (prevChrom == "") prevChrom = tmp(0)

			if (tmp(0) == prevChrom){
				gapList = (tmp(1).toInt - 1000,tmp(2).toInt + 1000) :: gapList
			} else {
				gapsChrom += prevChrom -> gapList.sortBy(_._1).toArray
				//println( prevChrom + " > " + tmp(0) + " " + gapList.size + " " + gapsChrom(prevChrom).size)
				prevChrom = tmp(0)
				gapList = Nil
				gapList = (tmp(1).toInt - 1000,tmp(2).toInt + 1000) :: gapList
			}
		}
		gapsChrom += prevChrom -> gapList.sortBy(_._1).toArray
		gap_in.close


		val curDepth = new HashMap[String,Double]
		while (depths.ready){
			val cline = depths.readLine.split("\t")
			curDepth += cline(0) -> cline(1).toDouble
		}
		depths.close
		var workingPop = new HashSet[String]

		var lines = new HashMap[String, PSEUDO] 

		while (dataSet.hasNext){

			val F = dataSet.next
			val fID = F.toString.split("/").last.split("\\.").apply(0)
			workingPop += fID
			val cEvents = new BufferedReader(new FileReader(F))

			while (cEvents.ready){
				val cline = cEvents.readLine
				//println(cline)
				val cEventLine = cline.split("\t")

				// Add depth based filtering if have >16x then must have several reads of support
				//if (curDepth(fID) >= 10.0 && (((cEventLine(1).toInt + cEventLine(2).toInt) >= curDepth(fID)/1.25 ) || ((cEventLine(3).toInt + cEventLine(4).toInt) >= 4) )){//&& (cEventLine(3).toInt + cEventLine(4).toInt) >= (curDepth(fID)/6) ){
				//Prob based filtering def bootStrap(reads: Int, evid: Int, geno: Int, avg: Int): Double = {
				//val difIPP = bootStrap((curDepth(fID)*15).toInt,(cEventLine(1).toInt + cEventLine(2).toInt),expReads,expIPP)
				//if ( difIPP <= 0.1){//&& (cEventLine(3).toInt + cEventLine(4).toInt) >= (curDepth(fID)/6) ){
				//	val difSR = bootStrap((curDepth(fID)*15).toInt,(cEventLine(3).toInt + cEventLine(4).toInt),expReads,expSR)
				System.err.println(fID)
				if (curDepth.contains(fID)){
				//	if ((difIPP * difSR) <= 0.0025 ){
							if (! lines.contains(cEventLine(0))){
									val coords = cEventLine(0).split(":").map(k => k.split("-")).flatten
									lines += cEventLine(0) -> new PSEUDO(coords(0),coords(1).toInt,coords(2).toInt)
								}
							//System.err.println(cEventLine(0))
							lines(cEventLine(0)).fwdIP = cEventLine(1).toInt + lines(cEventLine(0)).fwdIP
							lines(cEventLine(0)).revIP = cEventLine(2).toInt + lines(cEventLine(0)).revIP
							lines(cEventLine(0)).fwdSR = cEventLine(3).toInt + lines(cEventLine(0)).fwdSR 
							lines(cEventLine(0)).revSR = cEventLine(4).toInt + lines(cEventLine(0)).revSR
							lines(cEventLine(0)).allBreaks = addBreaksInfo(lines(cEventLine(0)).allBreaks,cEventLine(5))
							lines(cEventLine(0)).difs = cEventLine(6).toInt :: lines(cEventLine(0)).difs
							lines(cEventLine(0)).addRMInfo(cEventLine(7).split(":").map(_.toInt))
							lines(cEventLine(0)).addEventInfo(cEventLine(8).split(":").map(_.split(",")))
							lines(cEventLine(0)).addGeno(fID,cEventLine(9))
							lines(cEventLine(0)).windowDepth = cEventLine(10).toInt :: lines(cEventLine(0)).windowDepth
							lines(cEventLine(0)).other = cEventLine(11).toInt :: lines(cEventLine(0)).other
							lines(cEventLine(0)).readsRatio = cEventLine(12).toDouble :: lines(cEventLine(0)).readsRatio
							lines(cEventLine(0)).etype = cEventLine(13) :: lines(cEventLine(0)).etype
							lines(cEventLine(0)).addBreaks(cEventLine(14))
							//lines(cEventLine(0)).addSize(cEventLine(15))
							lines(cEventLine(0)).carriers = fID :: lines(cEventLine(0)).carriers
							lines(cEventLine(0)).addMSR(cEventLine(17))
						//}				
					} else {/*
						if (curDepth(fID) < 10 && (cEventLine(1).toInt + cEventLine(2).toInt) >= (curDepth(fID) * 2) ){
							if (! lines.contains(cEventLine(0))){
								val coords = cEventLine(0).split(":").map(k => k.split("-")).flatten
								lines += cEventLine(0) -> new LINE(coords(0),coords(1).toInt,coords(2).toInt)
							}
									
							lines(cEventLine(0)).fwdIP = cEventLine(1).toInt + lines(cEventLine(0)).fwdIP
							lines(cEventLine(0)).revIP = cEventLine(2).toInt + lines(cEventLine(0)).revIP
							lines(cEventLine(0)).fwdSR = cEventLine(3).toInt + lines(cEventLine(0)).fwdSR 
							lines(cEventLine(0)).revSR = cEventLine(4).toInt + lines(cEventLine(0)).revSR
							lines(cEventLine(0)).allBreaks = addBreaksInfo(lines(cEventLine(0)).allBreaks,cEventLine(5))
							lines(cEventLine(0)).difs = cEventLine(6).toInt :: lines(cEventLine(0)).difs
							lines(cEventLine(0)).addRMInfo(cEventLine(7).split(":").map(_.toInt))
							lines(cEventLine(0)).addEventInfo(cEventLine(8).split(":").map(_.split(",")))
							lines(cEventLine(0)).addGeno(fID,cEventLine(9))
							lines(cEventLine(0)).windowDepth = cEventLine(10).toInt :: lines(cEventLine(0)).windowDepth
							lines(cEventLine(0)).other = cEventLine(11).toInt :: lines(cEventLine(0)).other
							lines(cEventLine(0)).readsRatio = cEventLine(12).toDouble :: lines(cEventLine(0)).readsRatio
							lines(cEventLine(0)).etype = cEventLine(13) :: lines(cEventLine(0)).etype
							lines(cEventLine(0)).addBreaks(cEventLine(14))
							lines(cEventLine(0)).addSize(cEventLine(15))
							lines(cEventLine(0)).carriers = fID :: lines(cEventLine(0)).carriers	

						}
						*/
					}


			}
		}//Load data


		def isGAP (pos: Int, curSearch: Array[Tuple2[Int,Int]]) : Tuple3[Boolean,String,String] = {
			if (pos > 1){
				val half = scala.math.floor(curSearch.length/2).toInt
				if (curSearch.length == 0){
						//println(pos + " no gap")
						(false,"","")
				} else {
					if (curSearch(half)._1 <= pos && curSearch(half)._2 >= pos) {
						//println(pos + " gap " + curSearch(half)._1  + " " + curSearch(half)._2)
						(true,"gap","gap")
						} else {
							if (curSearch.length <= 1) {
								//println(pos + " nogap " + curSearch(half)._1  + " " + curSearch(half)._2)
								(false,"","")
							} else {
							if (curSearch(half)._1 > pos) isGAP(pos,curSearch.slice(0,half)) else isGAP(pos,curSearch.slice(half + 1,curSearch.length))
						}
					}
				}
				} else {
					(false,"","")
				}
		}
		var gtErrors = 0

		var newlines = new HashMap[String, PSEUDO] 
		/* Proper merger */
		for (elm <- lines.keys){
			for (elmT <- lines.keys){
				if (elm != elmT && lines(elm).chrom == lines(elmT).chrom && lines(elm).allBreaks.keys.toArray.map(a => lines(elmT).allBreaks.keys.toArray.contains(a)).contains(true)){
					//Same so Merge
							lines(elm).fwdIP += lines(elmT).fwdIP
							lines(elm).revIP += lines(elmT).revIP
							lines(elm).fwdSR += lines(elmT).fwdSR 
							lines(elm).revSR += lines(elmT).revSR

							lines(elm).difs = lines(elm).difs ++ lines(elmT).difs
							lines(elm).windowDepth = lines(elm).windowDepth ++ lines(elmT).windowDepth
							lines(elm).other = lines(elm).other ++ lines(elmT).other
							lines(elm).etype = lines(elm).etype ++ lines(elmT).etype
							lines(elm).carriers = lines(elm).carriers ++ lines(elmT).carriers
							lines(elm).readsRatio = lines(elm).readsRatio ++ lines(elmT).readsRatio

							//HashMap
							//lines(cEventLine(0)).allBreaks = addBreaksInfo(lines(cEventLine(0)).allBreaks,cEventLine(5))
							lines(elm).allBreaks = lines(elm).allBreaks ++ lines(elmT).allBreaks.map{ case (k,v) => k -> (v + lines(elm).allBreaks.getOrElse(k,0)) }

							//lines(cEventLine(0)).addRMInfo(cEventLine(7).split(":").map(_.toInt))
							lines(elm).rPmP += lines(elmT).rPmP
							lines(elm).rPmN += lines(elmT).rPmN
							lines(elm).rNmP += lines(elmT).rNmP
							lines(elm).rNmN += lines(elmT).rNmN


							//lines(cEventLine(0)).addEventInfo(cEventLine(8).split(":").map(_.split(",")))
							lines(elm).events = lines(elm).mergeData(lines(elm).events,lines(elmT).events)

							//lines(cEventLine(0)).addGeno(fID,cEventLine(9))
							lines(elm).geno = lines(elm).mergeData(lines(elm).geno,lines(elmT).geno)
							
							//lines(cEventLine(0)).addBreaks(cEventLine(14))
							lines(elm).best5p = lines(elm).best5p  ++ lines(elmT).best5p 
							lines(elm).best3p = lines(elm).best3p  ++ lines(elmT).best3p


							//lines(cEventLine(0)).addMSR(cEventLine(17))
							lines(elm).mateSR += lines(elmT).mateSR

						//Delete merged
					lines.remove(elmT)
				}

			}
		}



		println(s"Pos\t5'_IPP\t3'_IPP\t5'_SR\t3'_SR\tMate_SR\tCarriers\tCarrier%\tR+M+\tR+M-\tR-M+\tR-M-\tDenovo\tL1_Size\tTrio\tMendelian\tBest_Breakpoints\tBest_Event\tEvent_support\tGap?\tGenes\t")
		for (elm <- lines.values){
			val bestBreaks = if (elm.allBreaks.size >= 2) elm.allBreaks.toList.sortBy(-_._2).slice(0,2).sortBy(_._1) else if (elm.allBreaks.size == 1)  elm.allBreaks.toList ::: List((0,0)) else List((0,0),(0,0))
			//Check if near Gap
			var cleanOr = if ( ((elm.rPmP + elm.rNmN + 1)/(elm.rPmN + elm.rNmP + 1)) > 1 || ((elm.rPmN + elm.rNmP + 1)/(elm.rPmP + elm.rNmN + 1)) > 1 ) true else false
			
/* Total reads Sig */
				var depth = 0
				for (i <- elm.carriers) depth += curDepth(i).toInt
				val totalSeq = 12 * depth
				val splitP = bootStrap(totalSeq,elm.fwdSR + elm.revSR,expReads,expSR)
				val ippP = bootStrap(totalSeq,elm.fwdIP + elm.revIP,expReads,expIPP)
				System.err.println(s"${elm.chrom}:${elm.start}-${elm.end} IPP ${elm.fwdIP},${elm.revIP} P=${ippP}\tSR ${elm.fwdSR},${elm.revSR} P=${splitP}\tRatioIP/SR ${((elm.fwdIP + elm.revIP)/(elm.fwdSR + elm.revSR + 1))}")
			//if (!isGAP(bestBreaks(0)._1,gapsChrom(elm.chrom))._1 && (elm.fwdSR + elm.revSR) >= 4 && (elm.fwdIP + elm.revIP) >= 10 && cleanOr && ((elm.fwdIP + elm.revIP)/(elm.fwdSR + elm.revSR)) >= 1) {
				//if (!isGAP(bestBreaks(0)._1,gapsChrom(elm.chrom))._1 && ippP <= 0.05 && splitP <= 0.05 && cleanOr && ((elm.fwdIP + elm.revIP)/(elm.fwdSR + elm.revSR + 1.0)) >= 1.0) {
				if (ippP <= 0.05 && splitP <= 0.05 && cleanOr && ((elm.fwdIP + elm.revIP)/(elm.fwdSR + elm.revSR + 1.0)) >= 1.0) {
			//if (!isGAP(bestBreaks(0)._1,gapsChrom(elm.chrom))._1 && (elm.fwdSR + elm.revSR) >= 4 && (elm.fwdIP + elm.revIP) >= 10 && cleanOr && ((elm.fwdIP + elm.revIP)/(elm.fwdSR + elm.revSR)) >= 1) {
				var denovos : List[String] = Nil
				var denoNum, mendNum = 0
				var trio = ""
				for (carrier <- elm.carriers){
					if (ped.contains(carrier)){
							if (workingPop.contains(carrier)){
									if (workingPop.contains(ped(carrier).apply(2)) && workingPop.contains(ped(carrier).apply(3)) && !elm.carriers.contains(ped(carrier).apply(2)) && !elm.carriers.contains(ped(carrier).apply(3)) ){
										denovos = carrier :: denovos
										denoNum += 1
										if (denoNum == 1) trio = carrier + "_" + ped(carrier).apply(2) + "_" + ped(carrier).apply(3)
									} else {
										if (workingPop.contains(ped(carrier).apply(2)) && workingPop.contains(ped(carrier).apply(3)) && (elm.carriers.contains(ped(carrier).apply(2)) || elm.carriers.contains(ped(carrier).apply(3)))) {
											mendNum += 1
										} else {
											if ( (workingPop.contains(ped(carrier).apply(2)) && elm.carriers.contains(ped(carrier).apply(2)) ) ||  ( workingPop.contains(ped(carrier).apply(3)) && elm.carriers.contains(ped(carrier).apply(3)))) {
												mendNum += 1
											}
										}
									}
								} else {
									System.err.println(s"Carrier: ${carrier} missing from working Individuals")
								}
						} else {
							System.err.println(s"Carrier: ${carrier} missing from Pedigree")
						}
				}
				if (denoNum > 1) {
					gtErrors += denoNum
					denoNum = 0

				}

				val bestClass = if (elm.events.size >= 1) elm.events.toList.sortBy(s => -(s._2._1 + s._2._2)).apply(0) else ("ERROR",(0,0))
				//rPmP, rPmN, rNmP, rNmN 
				print(s"${elm.chrom}:${elm.start}-${elm.end}\t${elm.fwdIP}\t${elm.revIP}\t${elm.fwdSR}\t${elm.revSR}\t${elm.mateSR}\t${elm.carriers.size}\t${elm.carriers.size/workingPop.size.toFloat}\t${elm.rPmP}\t${elm.rPmN}\t${elm.rNmP}\t${elm.rNmN}\t${denoNum}\t\t${if (denoNum == 1) trio else "_"}\t${mendNum}\t${elm.chrom}:${bestBreaks(0)._1}-${bestBreaks(1)._1}\t")
				print(s"${bestClass._1}\t${bestClass._2._1}|${bestClass._2._2}\t${isGAP(bestBreaks(0)._1,gapsChrom(elm.chrom))._1}\t${isGene(elm.start,elm.chrom) + "_" + isGene(elm.end,elm.chrom)}")
				elm.carriers.foreach(s => print("\t" + s))
				print("\n")
			} else {
				System.err.println(s"${elm.chrom}:${elm.start}-${elm.end}\t${elm.fwdIP}\t${elm.revIP}\t${elm.fwdSR}\t${elm.revSR}\t${elm.rPmP}\t${elm.rPmN}\t${elm.rNmP}\t${elm.rNmN}\t${elm.chrom}:${bestBreaks(0)._1}-${bestBreaks(1)._1}\t${scala.math.abs(bestBreaks(1)._1 - bestBreaks(0)._1)}")
				System.err.println(s"GAP? ${!isGAP(bestBreaks(0)._1,gapsChrom(elm.chrom))._1}\t SR:${elm.fwdSR + elm.revSR}\tIPP:${elm.fwdIP + elm.revIP} less50:${scala.math.abs(bestBreaks(1)._1 - bestBreaks(0)._1) < 50}\tCleanOR:${cleanOr}")
			}
		}
		System.err.println(gtErrors + " rate est " + gtErrors/(130*lines.size).toFloat)
	}
}