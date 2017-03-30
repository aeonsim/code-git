import htsjdk.variant.vcf.VCFFileReader
import scala.collection.mutable.HashMap
import java.io._

System.err.println("scala -cp htsjdk.jar app.scala <inVCF> <outVCF> <ped> <probands> <maxAF> <minAF>")

val vcf = new VCFFileReader(new File(args(0)))
val builder = new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
val writer = builder.setOutputFile(args(1)).build
val ped = new BufferedReader(new FileReader(new File(args(2))))
val prob_in = new BufferedReader(new FileReader(new File(args(3))))

val maxAF = if (args.size >= 5) args(4).toDouble else 1.0
val minAF = if (args.size >= 6) args(5).toDouble else 0.0

var probands : List[String] = Nil

while (prob_in.ready){
	probands = prob_in.readLine :: probands
}

prob_in.close

var trios = new HashMap[String, Tuple2[String,String]]

while(ped.ready){
	val fam = ped.readLine.split("\t")
	if (probands.contains(fam(1))) trios += fam(1) -> Tuple2(fam(2),fam(3))
}

ped.close

writer.writeHeader(vcf.getFileHeader)

val vcfIT = vcf.iterator

while (vcfIT.hasNext){
	var keep = false
	val record = vcfIT.next
	val AF = if (record.getAttributeAsString("AF","").contains(",")) 99.0 else record.getAttributeAsDouble("AF",0.0)
	if (record.getAttributeAsDouble("MQ",0.0) == 60 && AF <= maxAF && AF >= minAF){

		var genos = record.getGenotypes
		
		for (fam <- trios){

			val pro = genos.get(fam._1)

			if (pro.isHet && pro.getDP >= 10){
				val pAR = pro.getAD.apply(1)/pro.getDP
				
				val sire = genos.get(fam._2._1)
				val dam = genos.get(fam._2._2)
				if(pAR >= 0.4 && (sire.getDP >= 10 && !sire.isHomRef) || (dam.getDP >= 10 && !dam.isHomRef)){
					keep = true
				}

			}
		}
	}

	if (keep){
		writer.add(record)
	}

}

writer.close





/*
val all = vcf.next
 val test = vcf.query("chr1",500000,1000000)

val mb = test.next
 mb.getAttributeAsInt("DP",0) >= 5000

 val gts = mb.getGenotypes

  val gt  = gts.get("10857209")

  gt.getDP
  gt.getGQ
  gt.getGenotypeString
  gt.isHomRef
  gt.isHomVar
  gt.isHet
  */