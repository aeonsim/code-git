import java.io._
import scala.collection.mutable.HashMap

val in = new BufferedReader(new FileReader(new File("Damona_Denovos_14th_Jan_2016.tab")))

val minAR = 0.3
val maxAR = 0.7
val minDPP = 0.2

var data = new HashMap[String, HashMap[String,Int]]

var tmp = in.readLine.split("\t")

while (in.ready){
	tmp = in.readLine.split("\t")
	if (tmp(24).toDouble >= minAR && tmp(24).toDouble <= maxAR && tmp(10).toDouble >= minDPP && tmp(4).size == 1 && tmp(5).size == 1){
		val event = tmp(4) + ">" + tmp(5)
		if (data.contains(tmp(7))){
			if (data(tmp(7)).contains(event)){
				data(tmp(7))(event) += 1
			} else {
				data(tmp(7)) += event -> 1
			}
		} else {
			data += tmp(7) -> HashMap(("C>A",0),("C>T",0),("C>G",0),("T>A",0),("T>C",0),("T>G",0),("G>A",0),("G>T",0),("G>C",0),("A>T",0),("A>C",0),("A>G",0))  
		}
	}

}

for (item <- data){
	val ts =  item._2("C>T") + item._2("G>A") + item._2("T>C") + item._2("A>G")
	val tv =  item._2("C>A") + item._2("A>C") + item._2("C>G") + item._2("G>C") + item._2("A>T") + item._2("T>A") + item._2("T>G") + item._2("G>T")
	println(item._1 + " ts/tv = " + (if (ts > 0 && tv > 0) ts/tv.toFloat) + " " + ts + " / " + tv)
}