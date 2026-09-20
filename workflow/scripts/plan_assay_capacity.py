#!/usr/bin/env python3
"""Calculate a prospective assay's well, image and instrument-capacity requirements.

These are planning arithmetic and explicit assumptions, not laboratory measurements,
power estimates, drug doses or evidence of biological efficacy.
"""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path


def positive_integer(value, label):
    if type(value) is not int or value < 1:
        raise ValueError(label + ' must be a positive integer')
    return value


def calculate(config):
    image = config['acquisition']
    names = ('duration_hours', 'cadence_minutes', 'fields_per_well', 'z_planes',
             'channels', 'width_pixels', 'height_pixels', 'stored_bytes_per_pixel')
    for name in names:
        positive_integer(image[name], name)
    duration_minutes = image['duration_hours'] * 60
    if duration_minutes % image['cadence_minutes']:
        raise ValueError('Duration must be exactly divisible by cadence for this inclusive-frame model')
    frames = duration_minutes // image['cadence_minutes'] + 1
    planes_per_well = image['fields_per_well'] * image['z_planes'] * image['channels']
    bytes_per_plane = image['width_pixels'] * image['height_pixels'] * image['stored_bytes_per_pixel']
    raw_bytes_per_well = frames * planes_per_well * bytes_per_plane
    copies = positive_integer(config['retained_raw_copies'], 'retained_raw_copies')
    capacities = config['simultaneous_well_capacities']
    if not capacities or len(set(capacities)) != len(capacities):
        raise ValueError('Supply distinct capacity scenarios')
    scenarios = []
    for capacity in capacities:
        positive_integer(capacity, 'simultaneous_well_capacity')
        planes = capacity * planes_per_well
        scenarios.append({'simultaneous_wells': capacity, 'image_planes_per_timepoint': planes,
                          'total_seconds_available_per_image_plane_including_all_overhead': image['cadence_minutes']*60/planes,
                          'raw_bytes_per_full_session': capacity*raw_bytes_per_well})
    stages = []
    for stage in config['stages']:
        factors = ('model_classes', 'founder_families', 'culture_preparations', 'conditions', 'technical_wells')
        for name in factors:
            positive_integer(stage[name], stage['name'] + '.' + name)
        total = math.prod(stage[name] for name in factors)
        per_preparation = total // stage['culture_preparations']
        instrument = []
        for capacity in capacities:
            sessions = stage['culture_preparations'] * math.ceil(per_preparation/capacity)
            instrument.append({'simultaneous_wells': capacity,
                               'minimum_sessions_preserving_separate_culture_preparations': sessions,
                               'minimum_acute_instrument_hours': sessions*image['duration_hours'],
                               'minimum_acute_instrument_days': sessions*image['duration_hours']/24})
        stages.append({'name':stage['name'], 'design':stage, 'wells':total,
                       'wells_per_culture_preparation':per_preparation,
                       'uncompressed_pixel_bytes':total*raw_bytes_per_well,
                       'raw_storage_bytes_with_copies':copies*total*raw_bytes_per_well,
                       'instrument_capacity':instrument})
    return {'kind':'prospective_planning_not_observed_experiment', 'configuration':config,
            'inclusive_frames_per_field_plane_channel':frames,
            'bytes_per_image_plane':bytes_per_plane, 'uncompressed_pixel_bytes_per_well':raw_bytes_per_well,
            'capacity_scenarios':scenarios, 'stages':stages,
            'limits':['No laboratory access, material availability, instrument throughput, efficacy or statistical power established.',
                      'Raw pixel counts exclude metadata, derived images, segmentation, temporary storage and spare capacity. No compression gain assumed.',
                      'Session counts are lower bounds, not a schedule. Setup, calibration, between-block effects, culture/model construction, washout and repeat failures require additional time.',
                      'Fields, frames and technical wells do not create independent biological replicates. The model assumes the specified founder families and culture preparations are actually available.',
                      'Three nonzero exposure conditions in the example are planning slots, not selected concentrations or a dosing recommendation.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Preserve previous plans; select a new output path')
    result = calculate(json.loads(args.config.read_text()))
    result['calculated_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['config_sha256'] = hashlib.sha256(args.config.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'output':str(args.output),'frames':result['inclusive_frames_per_field_plane_channel'],
                      'stages':[{k:s[k] for k in ('name','wells','uncompressed_pixel_bytes','raw_storage_bytes_with_copies')} for s in result['stages']]}))


if __name__ == '__main__':
    main()
