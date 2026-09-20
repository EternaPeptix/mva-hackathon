// Render a local research-draft pitch with native macOS media frameworks.
// Source images and speech timing are verified separately; no network or publication.
import Foundation
import AVFoundation
import AppKit
import CoreGraphics
import CoreVideo
import ImageIO

var createdOutput: URL? = nil
struct Segment: Decodable { let start_frame: Int; let end_frame: Int; let image: String }
struct Scene: Decodable { let id: String; let start_seconds: Double; let end_seconds: Double }
struct Plan: Decodable {
    let width: Int; let height: Int; let fps: Int; let duration_seconds: Int; let frames: Int
    let draft: Bool; let audio: String; let segments: [Segment]; let scenes: [Scene]
}
enum Failure: Error { case invalid(String) }
func require(_ condition: Bool, _ message: String) throws { if !condition { throw Failure.invalid(message) } }
func writeJSON(_ value: Any, _ url: URL) throws {
    let data = try JSONSerialization.data(withJSONObject: value, options: [.prettyPrinted, .sortedKeys])
    try data.write(to: url, options: .atomic)
}
func image(_ path: String) throws -> CGImage {
    guard let source = CGImageSourceCreateWithURL(URL(fileURLWithPath: path) as CFURL, nil),
          let image = CGImageSourceCreateImageAtIndex(source, 0, nil) else { throw Failure.invalid("Unreadable image") }
    return image
}
func buffer(_ image: CGImage, _ width: Int, _ height: Int) throws -> CVPixelBuffer {
    try require(image.width == width && image.height == height, "Image dimensions differ")
    var pixel: CVPixelBuffer?
    let attrs: [String: Any] = [kCVPixelBufferCGImageCompatibilityKey as String: true,
                               kCVPixelBufferCGBitmapContextCompatibilityKey as String: true]
    let result = CVPixelBufferCreate(kCFAllocatorDefault, width, height, kCVPixelFormatType_32ARGB, attrs as CFDictionary, &pixel)
    try require(result == kCVReturnSuccess && pixel != nil, "Pixel buffer creation failed")
    let p = pixel!; CVPixelBufferLockBaseAddress(p, [])
    defer { CVPixelBufferUnlockBaseAddress(p, []) }
    guard let context = CGContext(data: CVPixelBufferGetBaseAddress(p), width: width, height: height,
                                  bitsPerComponent: 8, bytesPerRow: CVPixelBufferGetBytesPerRow(p),
                                  space: CGColorSpaceCreateDeviceRGB(),
                                  bitmapInfo: CGImageAlphaInfo.noneSkipFirst.rawValue | CGBitmapInfo.byteOrder32Big.rawValue)
    else { throw Failure.invalid("Pixel context creation failed") }
    context.draw(image, in: CGRect(x: 0, y: 0, width: width, height: height))
    return p
}
func png(_ cg: CGImage, _ url: URL) throws {
    guard let destination = CGImageDestinationCreateWithURL(url as CFURL, "public.png" as CFString, 1, nil)
    else { throw Failure.invalid("PNG destination failed") }
    CGImageDestinationAddImage(destination, cg, nil)
    try require(CGImageDestinationFinalize(destination), "PNG write failed")
}

func run(planURL: URL, output: URL) throws {
    let start = Date(); let fm = FileManager.default
    try require(!fm.fileExists(atPath: output.path), "Preserve existing output; choose a new render directory")
    try fm.createDirectory(at: output.deletingLastPathComponent(), withIntermediateDirectories: true)
    try fm.createDirectory(at: output, withIntermediateDirectories: false)
    createdOutput = output
    let plan = try JSONDecoder().decode(Plan.self, from: Data(contentsOf: planURL))
    try require(plan.draft && plan.duration_seconds == 180 && plan.width == 1920 && plan.height == 1080 && plan.fps == 24 && plan.frames == 4320, "Unexpected draft specification")
    try require(plan.segments.first?.start_frame == 0 && plan.segments.last?.end_frame == plan.frames, "Incomplete frame plan")
    for (a,b) in zip(plan.segments, plan.segments.dropFirst()) { try require(a.end_frame == b.start_frame, "Frame gap or overlap") }
    try writeJSON(["status":"running", "stage":"encode original schematic frames", "started_at":ISO8601DateFormatter().string(from:start)], output.appendingPathComponent("status.json"))
    let silent = output.appendingPathComponent("silent-video.mp4")
    let writer = try AVAssetWriter(outputURL: silent, fileType: .mp4)
    let settings: [String:Any] = [AVVideoCodecKey: AVVideoCodecType.h264, AVVideoWidthKey: plan.width, AVVideoHeightKey: plan.height,
                                AVVideoCompressionPropertiesKey:[AVVideoAverageBitRateKey: 1800000, AVVideoExpectedSourceFrameRateKey: plan.fps,
                                                                 AVVideoMaxKeyFrameIntervalKey: 48, AVVideoProfileLevelKey: AVVideoProfileLevelH264HighAutoLevel]]
    let input = AVAssetWriterInput(mediaType: .video, outputSettings: settings)
    input.expectsMediaDataInRealTime = false
    let adaptor = AVAssetWriterInputPixelBufferAdaptor(assetWriterInput: input, sourcePixelBufferAttributes: nil)
    try require(writer.canAdd(input), "Cannot add video input"); writer.add(input)
    try require(writer.startWriting(), "Cannot start writing: \(String(describing: writer.error))")
    writer.startSession(atSourceTime: .zero)
    var buffers: [String:CVPixelBuffer] = [:]
    for segment in plan.segments {
        if buffers[segment.image] == nil { buffers[segment.image] = try buffer(image(segment.image), plan.width, plan.height) }
        for frame in segment.start_frame..<segment.end_frame {
            while !input.isReadyForMoreMediaData {
                try require(writer.status == .writing, "Video encoder stopped: \(String(describing: writer.error))")
                Thread.sleep(forTimeInterval: 0.003)
            }
            try require(adaptor.append(buffers[segment.image]!, withPresentationTime: CMTime(value: Int64(frame), timescale: Int32(plan.fps))), "Frame append failed")
        }
    }
    let duration = CMTime(value: 180, timescale: 1)
    writer.endSession(atSourceTime: duration); input.markAsFinished()
    let finish = DispatchSemaphore(value: 0); writer.finishWriting { finish.signal() }; finish.wait()
    try require(writer.status == .completed, "Video encoder did not finish: \(String(describing:writer.error))")
    buffers.removeAll()
    try writeJSON(["status":"running", "stage":"combine narration and video"], output.appendingPathComponent("status.json"))
    let video = AVURLAsset(url: silent); let audio = AVURLAsset(url: URL(fileURLWithPath: plan.audio))
    let composition = AVMutableComposition()
    guard let vsrc = video.tracks(withMediaType: .video).first, let asrc = audio.tracks(withMediaType: .audio).first,
          let vdst = composition.addMutableTrack(withMediaType: .video, preferredTrackID: kCMPersistentTrackID_Invalid),
          let adst = composition.addMutableTrack(withMediaType: .audio, preferredTrackID: kCMPersistentTrackID_Invalid)
    else { throw Failure.invalid("Missing media tracks") }
    try require(abs(audio.duration.seconds - 180) < 0.001, "Narration duration changed")
    let range = CMTimeRange(start: .zero, duration: duration)
    try vdst.insertTimeRange(range, of: vsrc, at: .zero)
    try adst.insertTimeRange(range, of: asrc, at: .zero)
    guard let exporter = AVAssetExportSession(asset: composition, presetName: AVAssetExportPreset1920x1080)
    else { throw Failure.invalid("Cannot create exporter") }
    let final = output.appendingPathComponent("jeweled_track2_pitch_draft_20260919.mp4")
    exporter.outputURL = final; exporter.outputFileType = .mp4; exporter.shouldOptimizeForNetworkUse = true; exporter.timeRange = range
    let exported = DispatchSemaphore(value:0); exporter.exportAsynchronously { exported.signal() }; exported.wait()
    try require(exporter.status == .completed, "Final export failed: \(String(describing:exporter.error))")
    try writeJSON(["status":"running", "stage":"decode and inspect the completed movie"], output.appendingPathComponent("status.json"))
    let asset = AVURLAsset(url: final)
    guard let track = asset.tracks(withMediaType:.video).first, let audioTrack = asset.tracks(withMediaType:.audio).first else { throw Failure.invalid("Exported tracks missing") }
    try require(abs(asset.duration.seconds - 180) < 0.001, "Final movie is not exactly 180 seconds")
    let format = track.formatDescriptions[0] as! CMFormatDescription
    let dimensions = CMVideoFormatDescriptionGetDimensions(format)
    try require(dimensions.width == 1920 && dimensions.height == 1080, "Final dimensions changed")
    let reader = try AVAssetReader(asset:asset); let stream = AVAssetReaderTrackOutput(track:track, outputSettings:nil)
    reader.add(stream); try require(reader.startReading(), "Cannot inspect encoded frames")
    var frames = 0; var buffersRead = 0; var emptyBuffers = 0; var firstPTS: Double?; var lastPTS: Double = -1
    while let sample = stream.copyNextSampleBuffer() {
        buffersRead += 1
        let count = CMSampleBufferGetNumSamples(sample)
        if count == 0 { emptyBuffers += 1; continue }
        let timestamp = CMSampleBufferGetPresentationTimeStamp(sample).seconds
        firstPTS = min(firstPTS ?? timestamp, timestamp); lastPTS = max(lastPTS,timestamp); frames += count
    }
    try require(reader.status == .completed && frames == 4320, "Frame count incomplete")
    // Decode the audio as signed little-endian PCM for an independent waveform comparison.
    let audioReader = try AVAssetReader(asset:asset)
    let audioOutput = AVAssetReaderTrackOutput(track:audioTrack, outputSettings:[AVFormatIDKey:kAudioFormatLinearPCM,AVLinearPCMBitDepthKey:16,AVLinearPCMIsFloatKey:false,AVLinearPCMIsBigEndianKey:false,AVLinearPCMIsNonInterleaved:false])
    audioReader.add(audioOutput); try require(audioReader.startReading(), "Cannot decode final audio")
    let pcm = output.appendingPathComponent("decoded-audio.pcm")
    fm.createFile(atPath:pcm.path,contents:nil); let handle = try FileHandle(forWritingTo:pcm); var audioSamples = 0
    while let sample = audioOutput.copyNextSampleBuffer() {
        guard let block = CMSampleBufferGetDataBuffer(sample) else { throw Failure.invalid("Missing decoded audio buffer") }
        let length = CMBlockBufferGetDataLength(block); var data = Data(count:length)
        let status = data.withUnsafeMutableBytes { bytes in CMBlockBufferCopyDataBytes(block, atOffset:0,dataLength:length,destination:bytes.baseAddress!) }
        try require(status == kCMBlockBufferNoErr, "Audio copy failed"); try handle.write(contentsOf:data); audioSamples += CMSampleBufferGetNumSamples(sample)
    }
    try handle.close(); try require(audioReader.status == .completed, "Audio decode incomplete")
    let generator = AVAssetImageGenerator(asset:asset); generator.appliesPreferredTrackTransform = true
    generator.requestedTimeToleranceBefore = .zero; generator.requestedTimeToleranceAfter = .zero
    var screenshots: [[String:Any]] = []
    for (i,scene) in plan.scenes.enumerated() {
        let seconds = scene.start_seconds + 2; let requested = CMTime(seconds:seconds,preferredTimescale:24); var actual = CMTime.zero
        let cg = try generator.copyCGImage(at:requested,actualTime:&actual)
        let file = output.appendingPathComponent(String(format:"decoded-scene-%02d.png",i+1));try png(cg,file)
        screenshots.append(["scene":scene.id,"requested_seconds":requested.seconds,"actual_seconds":actual.seconds,"path":file.path])
    }
    let receipt: [String:Any] = ["status":"complete","completed_at":ISO8601DateFormatter().string(from:Date()),"elapsed_seconds":Date().timeIntervalSince(start),
                               "video":final.path,"duration_seconds":asset.duration.seconds,"width":dimensions.width,"height":dimensions.height,
                               "video_frames":frames,"compressed_sample_buffers":buffersRead,"empty_control_buffers":emptyBuffers,"nominal_frame_rate":track.nominalFrameRate,"first_video_PTS":firstPTS!,"last_video_PTS":lastPTS,
                               "audio_track_duration_seconds":audioTrack.timeRange.duration.seconds,"decoded_audio_samples":audioSamples,"decoded_audio_PCM":pcm.path,
                               "video_codec_fourCC":CMFormatDescriptionGetMediaSubType(format),"decoded_scenes":screenshots,
                               "scope":"Local narrated research draft, not a published submission or clinical result"]
    try writeJSON(receipt,output.appendingPathComponent("render.json"))
    try writeJSON(["status":"complete","completed_at":ISO8601DateFormatter().string(from:Date())],output.appendingPathComponent("status.json"))
    print(String(data:try JSONSerialization.data(withJSONObject:receipt,options:[.sortedKeys]),encoding:.utf8)!)
}

do {
    try require(CommandLine.arguments.count == 3,"Usage: render_track2_pitch PLAN_JSON NEW_OUTPUT_DIRECTORY")
    try run(planURL:URL(fileURLWithPath:CommandLine.arguments[1]),output:URL(fileURLWithPath:CommandLine.arguments[2]))
} catch {
    fputs("Render failed; existing artifacts retained: \(error)\n",stderr)
    if let out = createdOutput {
        try? writeJSON(["status":"failed","error":String(describing:error)],out.appendingPathComponent("status.json"))
    }
    exit(1)
}
