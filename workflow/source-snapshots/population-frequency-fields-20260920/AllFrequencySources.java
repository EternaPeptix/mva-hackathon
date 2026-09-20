// Read-only exact-key public resource audit. No variant normalization or alignment access.
import java.io.*;
import java.nio.file.*;
import java.util.*;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.protobuf.*;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.lirical.exomiser_db_adapter.MvStoreUtil;

class AllFrequencySources {
    static final ObjectMapper JSON = new ObjectMapper();
    static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }
    static String b64(byte[] bytes) { return Base64.getEncoder().encodeToString(bytes); }
    static boolean included(String name) { return name.startsWith("GNOMAD_E_") || name.startsWith("GNOMAD_G_"); }
    static Object scalar(Object value) {
        if (value instanceof Descriptors.EnumValueDescriptor e)
            return Map.of("name", e.getName(), "number", e.getNumber());
        if (value instanceof Float f && !Float.isFinite(f)) return f.toString();
        if (value instanceof Double d && !Double.isFinite(d)) return d.toString();
        if (value instanceof ByteString b) return Map.of("base64", b64(b.toByteArray()));
        return value;
    }
    static Map<String,Object> schema() {
        var descriptor = AlleleProto.Frequency.getDescriptor();
        var fields = new ArrayList<Map<String,Object>>();
        for (var field : descriptor.getFields()) {
            var item = new LinkedHashMap<String,Object>();
            item.put("name", field.getName()); item.put("number", field.getNumber());
            item.put("type", field.getType().toString()); item.put("repeated", field.isRepeated());
            item.put("has_presence", field.hasPresence());
            item.put("default_getter_value", scalar(field.getDefaultValue()));
            if (field.getJavaType() == Descriptors.FieldDescriptor.JavaType.ENUM)
                item.put("enum_values", field.getEnumType().getValues().stream().map(v -> Map.of("name",v.getName(),"number",v.getNumber())).toList());
            fields.add(item);
        }
        return Map.of("message_name", descriptor.getFullName(), "fields", fields,
                      "file_descriptor_proto_base64", b64(descriptor.getFile().toProto().toByteArray()),
                      "presence_limitation", "For implicit-presence scalar fields, zero/default and absent cannot be distinguished after protobuf decoding. Positive AN does not recover AC presence: AC getter zero may be omitted/default rather than explicitly recorded zero. Reflective nondefault membership is not a presence indicator. Frequency scalar units are not assumed; no AF is calculated from that field.");
    }
    static Map<String,Object> frequency(AlleleProto.Frequency f, int ordinal) {
        var row = new LinkedHashMap<String,Object>();
        String source = f.getFrequencySource().name();
        row.put("ordinal_1based", ordinal); row.put("source_name", source);
        row.put("source_number", f.getFrequencySourceValue()); row.put("included_by_current_prefix_policy", included(source));
        row.put("ac_getter", f.getAc()); row.put("an_getter", f.getAn()); row.put("hom_getter", f.getHom());
        row.put("frequency_scalar_getter_uninterpreted", scalar(f.getFrequency()));
        row.put("frequency_scalar_ieee754_hex", String.format("%08x", Float.floatToRawIntBits(f.getFrequency())));
        row.put("decoded_message_reserialized_base64", b64(f.toByteArray()));
        row.put("unknown_fields_reserialized_base64", b64(f.getUnknownFields().toByteArray()));
        row.put("wire_limitation", "Serialized bytes are the decoded message reserialized, not the original MVStore wire bytes; explicit default encodings may already have been lost.");
        var fields = new LinkedHashMap<String,Object>();
        for (var field : f.getDescriptorForType().getFields()) {
            var v = new LinkedHashMap<String,Object>();
            v.put("getter_value", scalar(f.getField(field))); v.put("has_presence", field.hasPresence());
            v.put("present", field.hasPresence() ? f.hasField(field) : null);
            v.put("in_reflective_nondefault_field_map", f.getAllFields().containsKey(field));
            fields.put(field.getName(), v);
        }
        row.put("fields", fields);
        int ac = f.getAc(), an = f.getAn();
        String state = ac < 0 || an < 0 || ac > an ? "invalid_count_getters" :
            an == 0 ? "zero_denominator_or_implicit_default_unresolved" :
            ac == 0 ? "zero_AC_with_positive_AN" : "positive_AC_with_positive_AN";
        row.put("count_getter_state", state);
        // Ratio is descriptive only and never derived from the protobuf frequency scalar.
        row.put("AC_over_AN_if_valid_positive_denominator", an > 0 && ac >= 0 && ac <= an ? (double)ac/an : null);
        row.put("hom_getter_negative", f.getHom() < 0);
        return row;
    }
    static Map<String,Object> record(String key, boolean found, List<AlleleProto.Frequency> list) {
        require(found || list.isEmpty(), "absent key cannot have frequency rows");
        var result = new LinkedHashMap<String,Object>();
        result.put("key", key); result.put("exact_key_found", found);
        var all = new ArrayList<Map<String,Object>>(); var projected = new ArrayList<Map<String,Object>>();
        for (int i=0; i<list.size(); i++) {
            var f=list.get(i); all.add(frequency(f,i+1));
            if (included(f.getFrequencySource().name()))
                projected.add(Map.of("source",f.getFrequencySource().name(),"ac",f.getAc(),"an",f.getAn(),"hom",f.getHom()));
        }
        result.put("all_stored_frequency_rows_in_order", all);
        result.put("current_probe_projection", Map.of("key",key,"allele_found",found,"frequencies",projected));
        result.put("stored_frequency_row_count", all.size()); result.put("selected_prefix_row_count", projected.size());
        result.put("excluded_prefix_row_count", all.size()-projected.size());
        result.put("empty_export_explanation", !found ? "exact_key_not_found_representation_or_resource_reason_unresolved" :
            list.isEmpty() ? "found_key_has_no_stored_frequency_rows" :
            projected.isEmpty() ? "all_stored_frequency_rows_excluded_by_current_prefix_policy" : "selected_prefix_rows_present");
        return result;
    }
    static List<String[]> queries(Path path) throws IOException {
        var result = new ArrayList<String[]>(); var seen = new HashSet<String>();
        for (String line : Files.readAllLines(path)) {
            String[] q=line.split("\t",-1);
            require(q.length==4, "query must have four tab-separated fields");
            require(q[0].matches("[1-9][0-9]?") && Integer.parseInt(q[0])<=22, "invalid autosome");
            require(q[1].matches("[1-9][0-9]*") && Long.parseLong(q[1])<=Integer.MAX_VALUE, "invalid position");
            require(q[2].matches("[ACGT]+") && q[3].matches("[ACGT]+") && !q[2].equals(q[3]), "invalid concrete replacement");
            require(q[2].length()<=1000 && q[3].length()<=1000, "unbounded allele");
            require(seen.add(String.join(":",q)), "duplicate query"); result.add(q);
        }
        require(result.size()==25, "exactly 25 fixed keys required"); return result;
    }
    static List<String> selfTest() throws Exception {
        var passed = new ArrayList<String>();
        var e=AlleleProto.Frequency.newBuilder().setFrequencySource(AlleleProto.FrequencySource.GNOMAD_E_AFR).setAc(0).setAn(100).build();
        require(included("GNOMAD_E_AFR") && included("GNOMAD_G_MID") && !included("KG") && !included("ALFA_AFR"), "prefix policy"); passed.add("included_and_excluded_source_names");
        require(frequency(e,1).get("count_getter_state").equals("zero_AC_with_positive_AN"), "getter zero"); passed.add("zero_AC_getter_positive_AN_preserves_presence_uncertainty");
        require(frequency(e,1).get("AC_over_AN_if_valid_positive_denominator").equals(0.0), "zero ratio"); passed.add("zero_ratio_not_null");
        var absent=AlleleProto.Frequency.newBuilder().build();
        var explicitZero=AlleleProto.Frequency.newBuilder().setAc(0).setAn(0).setHom(0).setFrequency(0).build();
        require(Arrays.equals(absent.toByteArray(),explicitZero.toByteArray()), "implicit default semantics changed"); passed.add("explicit_zero_and_absent_scalar_wire_indistinguishable");
        require(AlleleProto.Frequency.getDescriptor().getFields().stream().noneMatch(Descriptors.FieldDescriptor::hasPresence), "unexpected schema presence"); passed.add("implicit_presence_schema_qualified");
        require(frequency(absent,1).get("AC_over_AN_if_valid_positive_denominator")==null, "default must not be zero AF"); passed.add("zero_denominator_default_is_unresolved");
        require(frequency(e.toBuilder().setAc(5).build(),1).get("AC_over_AN_if_valid_positive_denominator").equals(.05), "positive ratio"); passed.add("positive_valid_counts_ratio");
        require(frequency(e.toBuilder().setAc(101).build(),1).get("count_getter_state").equals("invalid_count_getters"), "invalid AC"); passed.add("AC_above_AN_refused_as_frequency");
        require(frequency(e.toBuilder().setAc(-1).build(),1).get("AC_over_AN_if_valid_positive_denominator")==null, "negative AC"); passed.add("negative_counts_not_frequency");
        var scalarOnly=absent.toBuilder().setFrequency(3.5f).build();
        require(frequency(scalarOnly,1).get("AC_over_AN_if_valid_positive_denominator")==null, "scalar frequency promoted"); passed.add("frequency_scalar_never_substitutes_counts");
        var unknown=e.toBuilder().setFrequencySourceValue(999).build();
        require(frequency(unknown,1).get("source_number").equals(999) && frequency(unknown,1).get("source_name").equals("UNRECOGNIZED") && !included(unknown.getFrequencySource().name()), "unknown enum"); passed.add("unknown_enum_number_preserved");
        var extras=e.toBuilder().setUnknownFields(UnknownFieldSet.newBuilder().addField(999,UnknownFieldSet.Field.newBuilder().addVarint(42).build()).build()).build();
        require(!frequency(extras,1).get("unknown_fields_reserialized_base64").equals(""), "unknown fields lost"); passed.add("unknown_fields_preserved");
        require(record("1:1:A:C",false,List.of()).get("empty_export_explanation").equals("exact_key_not_found_representation_or_resource_reason_unresolved"), "absent key"); passed.add("absent_key_distinct_from_found_empty");
        require(record("1:1:A:C",true,List.of()).get("empty_export_explanation").equals("found_key_has_no_stored_frequency_rows"), "found empty"); passed.add("found_key_without_rows");
        var kg=e.toBuilder().setFrequencySource(AlleleProto.FrequencySource.KG).build();
        require(record("1:1:A:C",true,List.of(kg)).get("empty_export_explanation").equals("all_stored_frequency_rows_excluded_by_current_prefix_policy"), "excluded evidence"); passed.add("excluded_only_export_explanation");
        var both=record("1:1:A:C",true,List.of(kg,e,e));
        require(both.get("stored_frequency_row_count").equals(3) && both.get("selected_prefix_row_count").equals(2), "ordered duplicates"); passed.add("duplicate_rows_preserved_no_deduplication");
        var restored=AlleleProto.Frequency.parseFrom(Base64.getDecoder().decode((String)frequency(extras,1).get("decoded_message_reserialized_base64")));
        require(restored.equals(extras), "round trip"); passed.add("full_decoded_frequency_message_round_trip");
        var nan=e.toBuilder().setFrequency(Float.NaN).build();
        require(frequency(nan,1).get("frequency_scalar_getter_uninterpreted").equals("NaN"), "nonfinite scalar"); passed.add("nonfinite_scalar_preserved_without_invalid_JSON");
        JSON.writeValueAsString(record("1:1:A:C",true,List.of(nan,unknown,extras)));
        return passed;
    }
    static Map<String,String> origins() {
        var out=new LinkedHashMap<String,String>();
        for (Class<?> c : List.of(AlleleProto.class,MvStoreUtil.class,MVStore.class,Message.class,ObjectMapper.class))
            out.put(c.getName(),c.getProtectionDomain().getCodeSource().getLocation().toString());
        return out;
    }
    public static void main(String[] args) throws Exception {
        if (args.length==1 && args[0].equals("--self-test")) {
            var passed=selfTest();
            System.out.println(JSON.writeValueAsString(Map.of("status","complete","tests_passed",passed.size(),"tests_failed",0,"tests_skipped",0,"tests",passed,"schema",schema(),"class_origins",origins(),"database_opened",false)));
            return;
        }
        require(args.length==6 && args[0].equals("--database") && args[2].equals("--queries") && args[4].equals("--output"), "use --self-test or explicit database/queries/new-output");
        Path database=Path.of(args[1]), queryFile=Path.of(args[3]), output=Path.of(args[5]);
        require(Files.isRegularFile(database) && !Files.exists(output), "missing database or existing output");
        var queries=queries(queryFile); var records=new ArrayList<Map<String,Object>>();
        try (var store=new MVStore.Builder().fileName(database.toString()).readOnly().open()) {
            var alleles=MvStoreUtil.openAlleleMVMap(store);
            for (String[] q : queries) {
                var key=AlleleProto.AlleleKey.newBuilder().setChr(Integer.parseInt(q[0])).setPosition(Integer.parseInt(q[1])).setRef(q[2]).setAlt(q[3]).build();
                var value=alleles.get(key);
                records.add(record(String.join(":",q),value!=null,value==null ? List.of() : value.getFrequenciesList()));
            }
        }
        var result=Map.of("status","complete","scope","read-only exact-key public resource audit; no population-frequency reassignment","records",records,"schema",schema(),"class_origins",origins());
        Files.writeString(output,JSON.writerWithDefaultPrettyPrinter().writeValueAsString(result)+"\n",StandardOpenOption.CREATE_NEW);
    }
}
