//package recombDenovo

object recombDenovo {

	def main(args: Array[String]) : Unit = {

			import java.io._
			import scala.collection.mutable.HashMap
			if (args.size == 0 || args.contains("-h")){
				println("recombDenovo <denovo list>\n<denovo list> = Text file with tab separated Sample Name \t Chr \t Position")
				System.exit(1)
			}
			var inRecomb, notRecomb = 0
					val files = new File(".").listFiles()
					val denovoIn = new BufferedReader(new FileReader(args(0)))
			var denovos = new HashMap[String, HashMap[String,List[Tuple3[Int,String,String]]]]

					while (denovoIn.ready){
						var cline = denovoIn.readLine.split("\t")
								if (denovos.contains(cline(0))) {
									if (denovos(cline(0)).contains(cline(1))) {
										denovos(cline(0))(cline(1)) = (cline(2).toInt,cline(3),cline(4)) :: denovos(cline(0))(cline(1))
									} else {
										denovos(cline(0)) += cline(1) -> List((cline(2).toInt,cline(3),cline(4)))
									}
								}else{
									denovos += cline(0) -> HashMap(cline(1) -> List((cline(2).toInt,cline(3),cline(4))))
								}
					}
			var prev, blockstart = Array("")
					for (cfile <- files){
						if(cfile.toString.contains(".bed") && cfile.length >= 1000){
							//details(0) Child name, details(1) parentName-origin.bed 
							val details = cfile.toString.split("/")(1).split("_")
							var cIn = new BufferedReader(new FileReader(cfile))
							//println(cfile)
							while (cIn.ready){
								var cline = cIn.readLine.split(" ")
										if (cline(4).toInt >= 10 && denovos.contains(details(0))){
											if(cline(0) != prev(0)){
											    blockstart = cline
												prev = cline
											} else {
												if (cline(3) != prev(3)){
													var start = prev(2).toInt
															var end = cline(1).toInt
															var cAn = denovos(details(0))
															if (cAn.contains(cline(0))) {
																for (dnM <- cAn(cline(0))){
																	if (dnM._1 >= start && dnM._1 <= end){
																		println(details.foreach(s => print(s + " ")) + " Recombination Denovo " + cline(0) + ":" + start + "-" + end + "\tPosition " + dnM._1 + " Score " + dnM._2 + " Window Size " + (end-start) + " Origin " + dnM._3 )
																		inRecomb += 1
																	} //If within start/end
																}//for denovos in chr
															}//if Chr exists
													blockstart = cline
													prev = cline
												} else { // if same origin
													prev = cline

												}//else same chr

											}//else same chr

										}//If > 10

							}

						}//if
					}
			println("Total De novo Mutations in Recomb " + inRecomb)
	}//Main
}//Ob