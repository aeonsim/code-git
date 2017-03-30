import net.sf.samtools.util._
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap
import java.io._

/*

val args = Array("11","11off_downsampled.vcf.gz","PR1denovo.txt","fam.txt","NL288458773")

*/


//println("app <numkids> <vcf> <denovos> <family.tab> familyID")

/* chrom	pos*/
val inDenovo = new BufferedReader(new FileReader(new File(args(2))))
/* Pro	kid1	kid2	...*/
val inFam = new BufferedReader(new FileReader(new File(args(3))))
val in_VCF = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(1)))))

var denovos = new HashSet[String]
var data = new HashMap[Int,List[String]]
var fam = new HashMap[String,List[Int]]

var found, ngtkids = 0
val numKids = args(0).toInt

while (inDenovo.ready){
val tmp = inDenovo.readLine.split("\t")
denovos += s"${tmp(0)}:${tmp(1)}"
}
inDenovo.close

var line = in_VCF.readLine.split("\t")

while (line(0)(1) == '#') line = in_VCF.readLine.split("\t")

//while (inFam.ready){
val tmp = inFam.readLine.split("\t")
var kids : List[Int] = Nil
fam += tmp(0) -> kids
for (i <- 1 to (tmp.size - 1)){
fam(tmp(0)) = line.indexOf(tmp(i)) :: fam(tmp(0))
}
//}


while (in_VCF.ready){
line = in_VCF.readLine.split("\t")

if(denovos.contains(s"${line(0)}:${line(1)}")){
found += 1
var gtkid = false
//System.err.println(s"${line(0)}:${line(1)}")
for (i <- fam){
for (k <- i._2){
val gtinfo = line(k).split(":")
if (gtinfo(0).size == 3) if (gtinfo(0)(0) == '1' || gtinfo(0)(2) == '1'){
gtkid = true
if (data.contains(k)){
data(k) = s"${line(0)}:${line(1)}" :: data(k)
} else {
data += k -> List(s"${line(0)}:${line(1)}")
}
}

}
}
if (gtkid) ngtkids += 1
}

}

//System.err.println("###########")

var results : List[Int] = Nil

for (j <- 0 to 100000){
var oldKids : List[Int] = Nil
var udnms = new HashSet[String] 
val ran = scala.util.Random
val kidsl = fam(args(4))
while (oldKids.size < numKids){
val nextK = ran.nextInt(kidsl.size)
if (oldKids.contains(nextK)){
} else {
oldKids = nextK :: oldKids
val dnms = data(kidsl(nextK))
for (i <- dnms) udnms += i
} 
}
//for (l <- denovos) if (!udnms.contains(l)) System.err.println(l)
results = udnms.size :: results
}

println("Found\t" + found + "\twith_Kids\t" + ngtkids + "\tAvg_DNMs\t" + results.sum/results.size.toFloat + 
	"\tKids\t" + numKids +  "\tsensitivity_(kids)_=\t" + (results.sum/results.size.toDouble/ngtkids).toFloat + 
	"\tOverall_senstivity_=\t" + (results.sum/results.size.toDouble/found).toFloat + "\tmin\t" + results.min + "\tmax\t" + results.max)