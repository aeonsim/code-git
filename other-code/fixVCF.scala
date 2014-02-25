import java.io._
import net.sf.samtools.util._

val input = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream("GATK-502-sorted.full.vcf.gz"))))
val out = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream("GATK-POP-CLEANED.VCF.GZ")))

var line = input.readLine

while (line(1) == '#') {
out.write(line)
line = input.readLine

}
val size = line.split("\t").size
out.write(line)
line = input.readLine

while (input.ready){
val curarr = line.split("\t")
if (curarr.size == size){
out.write(line)
} else {
println("Size should be: " + size + "\t Current Line is: " + curarr.size + " Position was: " + curarr(0) + " " + curarr(1))
}

line = input.readLine

}