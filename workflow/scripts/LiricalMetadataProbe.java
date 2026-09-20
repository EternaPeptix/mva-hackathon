// Read-only query of public annotation databases with LIRICAL's bundled reader.
// Usage: java -cp lirical-cli.jar LiricalMetadataProbe.java variants.mv.db clinvar.mv.db
import java.util.*;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.lirical.exomiser_db_adapter.MvStoreUtil;
import com.fasterxml.jackson.databind.ObjectMapper;

class LiricalMetadataProbe {
    public static void main(String[] args) throws Exception {
        if (args.length != 2) throw new IllegalArgumentException("Expected allele and ClinVar DB paths");
        // Candidate coordinates, not an upload or iteration over a subject's genome.
        String[][] queries = {
            {"15", "40209701", "T", "G"}, {"15", "40220612", "T", "G"},
            {"6", "98917691", "T", "C"}, {"12", "47906043", "T", "C"},
            {"22", "17600027", "C", "T"}, {"22", "17600050", "T", "G"}
        };
        var output = new ArrayList<Map<String,Object>>();
        try (var alleles = new MVStore.Builder().fileName(args[0]).readOnly().open();
             var clinvar = new MVStore.Builder().fileName(args[1]).readOnly().open()) {
            var alleleMap = MvStoreUtil.openAlleleMVMap(alleles);
            var clinvarMap = MvStoreUtil.openClinVarMVMap(clinvar);
            for (var query : queries) {
                var key = AlleleProto.AlleleKey.newBuilder()
                    .setChr(Integer.parseInt(query[0])).setPosition(Integer.parseInt(query[1]))
                    .setRef(query[2]).setAlt(query[3]).build();
                var item = new LinkedHashMap<String,Object>();
                item.put("grch38_allele", String.join(":", query));
                var annotation = alleleMap.get(key);
                var cv = clinvarMap.get(key);
                item.put("allele_found", annotation != null);
                if (annotation != null) {
                    item.put("rs_id", annotation.getRsId());
                    item.put("pathogenicity_scores", annotation.getPathogenicityScoresList().stream()
                        .map(p -> Map.of("source", p.getPathogenicitySource().toString(),
                                        "score", p.getScore())).toList());
                    item.put("frequency_sources", annotation.getFrequenciesList().stream()
                        .map(f -> Map.of("source", f.getFrequencySource().toString(),
                            "af_field", f.getFrequency(), "ac", f.getAc(), "an", f.getAn(),
                            "hom", f.getHom())).toList());
                }
                item.put("clinvar_found", cv != null);
                if (cv != null) {
                    item.put("clinvar_variation_id", cv.getVariationId());
                    item.put("clinvar_primary_interpretation", cv.getPrimaryInterpretation().toString());
                    item.put("clinvar_review_status", cv.getReviewStatus().toString());
                }
                output.add(item);
            }
        }
        System.out.println(new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString(output));
    }
}
