// Decode one midpoint frame from every planned static segment of a local movie.
// This checks all planned images, but does not claim exhaustive frame-by-frame QA.
import Foundation
import AVFoundation
import ImageIO
import CoreGraphics

struct Segment: Decodable { let start_frame: Int; let end_frame: Int; let image: String }
struct Plan: Decodable { let fps: Int; let frames: Int; let segments: [Segment] }
enum InspectionFailure: Error { case invalid(String) }
func require(_ value: Bool, _ message: String) throws {
    if !value { throw InspectionFailure.invalid(message) }
}
func run(planURL: URL, movieURL: URL, output: URL) throws {
    let fm = FileManager.default
    try require(!fm.fileExists(atPath: output.path), "Existing inspection directory preserved")
    let plan = try JSONDecoder().decode(Plan.self, from: Data(contentsOf: planURL))
    try require(plan.fps > 0 && plan.frames > 0 && !plan.segments.isEmpty, "Empty frame plan")
    try require(plan.segments.first!.start_frame == 0 && plan.segments.last!.end_frame == plan.frames, "Incomplete endpoints")
    for (a,b) in zip(plan.segments, plan.segments.dropFirst()) {
        try require(a.end_frame == b.start_frame, "Segment gap or overlap")
    }
    for segment in plan.segments {
        try require(segment.start_frame >= 0 && segment.end_frame > segment.start_frame && segment.end_frame <= plan.frames, "Invalid segment")
    }
    let asset = AVURLAsset(url: movieURL)
    try require(abs(asset.duration.seconds - Double(plan.frames)/Double(plan.fps)) < 0.001, "Movie duration disagrees with plan")
    let generator = AVAssetImageGenerator(asset: asset)
    generator.appliesPreferredTrackTransform = true
    generator.requestedTimeToleranceBefore = .zero
    generator.requestedTimeToleranceAfter = .zero
    try fm.createDirectory(at: output, withIntermediateDirectories: false)
    var rows: [[String:Any]] = []
    for (index, segment) in plan.segments.enumerated() {
        let frame = (segment.start_frame + segment.end_frame - 1)/2
        let requested = CMTime(value: Int64(frame), timescale: Int32(plan.fps))
        var actual = CMTime.zero
        let image = try generator.copyCGImage(at: requested, actualTime: &actual)
        try require(abs(actual.seconds - requested.seconds) < 0.5/Double(plan.fps), "Decoded frame time changed")
        let path = output.appendingPathComponent(String(format: "segment-%02d.png", index+1))
        guard let target = CGImageDestinationCreateWithURL(path as CFURL, "public.png" as CFString, 1, nil) else {
            throw InspectionFailure.invalid("PNG output unavailable")
        }
        CGImageDestinationAddImage(target, image, nil)
        try require(CGImageDestinationFinalize(target), "PNG output incomplete")
        rows.append(["segment_index":index, "start_frame":segment.start_frame, "end_frame":segment.end_frame,
                     "sampled_frame":frame, "requested_seconds":requested.seconds, "actual_seconds":actual.seconds,
                     "source_image":segment.image, "decoded_image":path.path, "width":image.width, "height":image.height])
    }
    let receipt: [String:Any] = ["status":"complete", "completed_at_utc":ISO8601DateFormatter().string(from:Date()),
                               "movie":movieURL.path, "segments":rows,
                               "scope":"One decoded midpoint per planned static segment; not every encoded frame"]
    let data = try JSONSerialization.data(withJSONObject:receipt,options:[.prettyPrinted,.sortedKeys])
    try data.write(to:output.appendingPathComponent("inspection.json"),options:.atomic)
    print("Decoded \(rows.count) planned segment midpoints")
}
do {
    try require(CommandLine.arguments.count == 4, "Usage: inspect_track2_pitch PLAN_JSON MOVIE NEW_OUTPUT_DIRECTORY")
    try run(planURL:URL(fileURLWithPath:CommandLine.arguments[1]), movieURL:URL(fileURLWithPath:CommandLine.arguments[2]), output:URL(fileURLWithPath:CommandLine.arguments[3]))
} catch {
    fputs("Inspection failed; existing files preserved: \(error)\n",stderr)
    exit(1)
}
