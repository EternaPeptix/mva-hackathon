import java.nio.file.*;
import java.util.*;
import genepi.mut.objects.BasePosition;
import genepi.mut.objects.VariantLine;
import com.google.gson.GsonBuilder;

/** Fictional rank-switch example, not a patient likelihood or sensitivity test. */
public class MutserveAlleleRankProbe {
    public static void main(String[] args) throws Exception {
        if(args.length!=1 || Files.exists(Paths.get(args[0]))) throw new IllegalArgumentException();
        Map<String,Object> output=new LinkedHashMap<>();
        for(boolean extra: new boolean[]{false,true}) {
            BasePosition b=new BasePosition();b.setId("fictional");b.setPos(1000);
            b.addtFor(extra?950:980);b.addgFor(20);b.addaFor(extra?30:0);
            b.addtRev(980);b.addgRev(20);
            for(int i=0;i<(extra?950:980);i++)b.addtForQ((byte)40);
            for(int i=0;i<20;i++){b.addgForQ((byte)40);b.addgRevQ((byte)40);}
            for(int i=0;i<(extra?30:0);i++)b.addaForQ((byte)40);
            for(int i=0;i<980;i++)b.addtRevQ((byte)40);
            VariantLine line=new VariantLine();line.setRef('T');line.parseLine(b,.01,null);
            Map<String,Object> row=new LinkedHashMap<>();
            row.put("FWD",Map.of("T",extra?950:980,"G",20,"A",extra?30:0));
            row.put("REV",Map.of("T",980,"G",20));row.put("eligible_minor_bases",line.getMinors());
            output.put(extra?"extra_forward_allele":"two_alleles",row);
        }
        Files.writeString(Paths.get(args[0]),new GsonBuilder().setPrettyPrinting().create().toJson(output)+"\n");
    }
}
