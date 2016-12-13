import java.io._

val ts = Array("A>G","G>A","C>T","T>C")
val tv = Array("C>G","G>C","A>T","T>A","A>C","C>A","G>T","T>G")
val obs = 100000
val lowlim = 1.1


//modlowAF.vcf
val in = new BufferedReader(new FileReader(new File(args(0))))
//val in = new BufferedReader(new FileReader(new File("modlowAF.vcf")))

var data : List[String] = Nil

var line = in.readLine

while (line(0) == '#') line = in.readLine

while (in.ready){
	val tmp = line.split("\t")
	if (tmp(3).size == 1 && tmp(4).size == 1) data = s"${tmp(3)}>${tmp(4)}" :: data

	line = in.readLine
}

var low = 0
val adata = data.toArray

for (i <- 1 to obs){
	var numTS, numTV = 0
	val rand = new scala.util.Random
	for (y <- 1 to 202){
		val sample = adata(rand.nextInt(adata.size))
		//print(sample)
		if (ts.contains(sample)) numTS += 1
		if (tv.contains(sample)) numTV += 1
	}
	val tstv = numTS.toFloat / numTV.toFloat
	//println(s"TS:${numTS} TV:${numTV} TSTV:${tstv}")
	if (tstv <= lowlim) low += 1

}

println (s"After ${obs} obs ${low} tstv's below ${lowlim} p = " + low/obs.toFloat)