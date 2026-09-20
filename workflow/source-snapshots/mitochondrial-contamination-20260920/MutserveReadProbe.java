import java.nio.file.*;
import java.util.*;
import genepi.mut.pileup.BamAnalyser;
import genepi.mut.objects.BasePosition;
import htsjdk.samtools.*;
import com.google.gson.GsonBuilder;

/** Invokes the release JAR's actual single-read counter on fictional reads. */
public class MutserveReadProbe {
    static int depth(BasePosition b) {
        if (b == null) return 0;
        return b.getaFor()+b.getaRev()+b.getcFor()+b.getcRev()+b.getgFor()+b.getgRev()+b.gettFor()+b.gettRev();
    }
    public static void main(String[] args) throws Exception {
        if (args.length != 2 || Files.exists(Paths.get(args[1]))) throw new IllegalArgumentException();
        String[] names = {"primary", "secondary", "supplementary", "qc_failed", "duplicate", "mq255", "mq29", "as29", "missing_as", "bq19", "length25", "overlapping_pair"};
        Map<String,Integer> result = new LinkedHashMap<>();
        for (String name : names) {
            BamAnalyser a = new BamAnalyser("fictional", args[0], 20, 30, 30, false, "mtdna");
            SAMFileHeader h = new SAMFileHeader(); h.addSequence(new SAMSequenceRecord("chrM",16569));
            int len = name.equals("length25") ? 25 : 100;
            SAMRecord r = new SAMRecord(h);
            r.setReadName(name); r.setReferenceName("chrM"); r.setAlignmentStart(1000);
            r.setReadString("A".repeat(len)); r.setBaseQualityString((name.equals("bq19") ? "4" : "I").repeat(len));
            r.setCigarString(len+"M"); r.setMappingQuality(name.equals("mq255") ? 255 : name.equals("mq29") ? 29 : 60);
            if (!name.equals("missing_as")) r.setAttribute("AS",name.equals("as29") ? 29 : 100);
            if (name.equals("secondary")) r.setFlags(256);
            if (name.equals("supplementary")) r.setFlags(2048);
            if (name.equals("qc_failed")) r.setFlags(512);
            if (name.equals("duplicate")) r.setFlags(1024);
            if (name.equals("overlapping_pair")) r.setFlags(65);
            a.analyseRead(r,true,false);
            if (name.equals("overlapping_pair")) {r.setFlags(129); a.analyseRead(r,true,false);}
            result.put(name,depth(a.getCounts().get(1000)));
        }
        Files.writeString(Paths.get(args[1]),new GsonBuilder().setPrettyPrinting().create().toJson(result)+"\n");
    }
}
