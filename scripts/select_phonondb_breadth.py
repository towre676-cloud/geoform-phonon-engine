"""Select one NIMS record per previously unattempted reduced composition."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from math import gcd
from pathlib import Path
import re


def normalized(counts):
    divisor = gcd(*counts.values())
    return tuple(sorted((element, count // divisor) for element, count in counts.items()))


def composition(formula):
    parts = re.findall(r'([A-Z][a-z]?)([0-9]*)', formula)
    if not parts or ''.join(element + count for element, count in parts) != formula:
        return None
    counts = Counter()
    for element, count in parts:
        counts[element] += int(count or '1')
    return normalized(counts) if all(counts.values()) else None


def explicit_profile_compositions(browse_index):
    """Use explicit compositions for priority only, without assigning identities."""
    column = browse_index['columns'].index('parsed')
    for row in browse_index['rows']:
        parsed = row[column]
        if not parsed or parsed.get('status') != 'EXPLICIT':
            continue
        try:
            parts = (item.split(':') for item in parsed['key'].split('|'))
            counts = {element: int(count) for element, count in parts}
        except (ValueError, KeyError):
            continue
        if counts and all(count > 0 for count in counts.values()):
            yield normalized(counts)


def tie_break(value):
    return hashlib.sha256(('GeoForm NIMS wave 2 ' + str(value)).encode()).hexdigest()


def select(rows, prior_ids, profile_compositions, limit):
    seen = {composition(row['formula']) for row in rows if row['record_id'] in prior_ids}
    candidates = defaultdict(list)
    for row in rows:
        key = composition(row['formula'])
        if key is not None and key not in seen:
            candidates[key].append(row['record_id'])
    ranked = sorted(candidates, key=lambda key: (key not in profile_compositions, tie_break(key)))
    return [min(candidates[key], key=tie_break) for key in ranked[:limit]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_index', type=Path)
    parser.add_argument('prior_ids', type=Path)
    parser.add_argument('browse_index', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--limit', type=int, default=750)
    args = parser.parse_args()
    if args.limit < 1:
        parser.error('selection limit must be positive')
    rows = json.loads(args.source_index.read_text())['rows']
    prior = set(json.loads(args.prior_ids.read_text()))
    profile_compositions = set(explicit_profile_compositions(json.loads(args.browse_index.read_text())))
    ids = select(rows, prior, profile_compositions, args.limit)
    args.output.write_text(json.dumps(ids, indent=2) + '\n')


if __name__ == '__main__':
    main()
