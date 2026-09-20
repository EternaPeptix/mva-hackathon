import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;
import importer.VcfImporter;
import vcf.Sample;
import genepi.haplocheck.steps.contamination.VariantSplitter;
import com.google.gson.GsonBuilder;
public class HaplocheckParserProbe {
    public static void main(String[] args) throws Exception {
        if (args.length != 2) throw new IllegalArgumentException("input VCF and new JSON output required");
        if (Files.exists(Paths.get(args[1]))) throw new IllegalArgumentException("existing output preserved");
        java.util.HashMap<String, Sample> samples = new VcfImporter().load(new File(args[0]), false);
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("samples", new TreeMap<String, Sample>(samples));
        result.put("profiles", new VariantSplitter().split(samples));
        Files.write(Paths.get(args[1]), (new GsonBuilder().setPrettyPrinting().create().toJson(result)+"\n").getBytes(StandardCharsets.UTF_8));
    }
}
