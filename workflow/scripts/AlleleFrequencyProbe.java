// Local read-only frequency lookup. Input: tab-separated chrom, pos, REF, ALT.
// Output: one JSON object per line; no phenotype, ClinVar or predictor scores.
import java.io.*;
import java.util.*;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.lirical.exomiser_db_adapter.MvStoreUtil;
import com.fasterxml.jackson.databind.ObjectMapper;

class AlleleFrequencyProbe {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("Expected public variants.mv.db path");
        var json = new ObjectMapper();
        try (var store = new MVStore.Builder().fileName(args[0]).readOnly().open();
             var input = new BufferedReader(new InputStreamReader(System.in));
             var output = new BufferedWriter(new OutputStreamWriter(System.out))) {
            var alleles = MvStoreUtil.openAlleleMVMap(store);
            for (String line; (line = input.readLine()) != null;) {
                var q = line.split("\t", -1);
                if (q.length != 4) throw new IllegalArgumentException("Expected four query fields");
                int chrom = Integer.parseInt(q[0]), pos = Integer.parseInt(q[1]);
                if (chrom < 1 || chrom > 22 || pos < 1 || !q[2].matches("[ACGT]+") || !q[3].matches("[ACGT]+"))
                    throw new IllegalArgumentException("Invalid autosomal allele query");
                var key = AlleleProto.AlleleKey.newBuilder().setChr(chrom).setPosition(pos)
                    .setRef(q[2]).setAlt(q[3]).build();
                var value = alleles.get(key);
                var record = new LinkedHashMap<String,Object>();
                record.put("key", String.join(":", q));
                record.put("allele_found", value != null);
                var frequencies = new ArrayList<Map<String,Object>>();
                if (value != null) {
                    for (var f : value.getFrequenciesList()) {
                        var name = f.getFrequencySource().toString();
                        if (name.startsWith("GNOMAD_E_") || name.startsWith("GNOMAD_G_"))
                            frequencies.add(Map.of("source", name, "ac", f.getAc(), "an", f.getAn(), "hom", f.getHom()));
                    }
                }
                record.put("frequencies", frequencies);
                output.write(json.writeValueAsString(record));
                output.newLine();
            }
        }
    }
}
