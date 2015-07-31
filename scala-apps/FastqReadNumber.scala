import net.sf.samtools.util._
import java.io._
import org.apache.commons.io.FileUtils._
import scala.collection.JavaConversions._

val files = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("gz"),true).iterator
var long = 0L

for (f <- files){
val path = f.toString.split("/")
val details = path(path.size - 1).split("_")
val cfile = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(f))))

while(cfile.ready){
val x = cfile.readLine
long += 1
}
println(path(path.size - 1) + "\t" + details(0) + "\t" + details(1) + "\t" + details(2) + "\t" + long/4)
}