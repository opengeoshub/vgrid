from collections import defaultdict

def get_default_stats_structure():
    return {
        'count': 0,
        'sum': [],
        'min': [],
        'max': [],
        'mean': [],
        'median': [],
        'std': [],
        'var': [],
        'range': [],
        'values': []  # For variety, minority, majority
    }

def append_stats_value(h3_bins, h3_id, props, stats, category_field, numeric_field=None):
    category_value = props.get(category_field, "all") if category_field else "all"
    if h3_id not in h3_bins:
        h3_bins[h3_id] = defaultdict(get_default_stats_structure)

    if stats == 'count':
        h3_bins[h3_id][category_value]['count'] += 1
    elif stats in ['minority', 'majority', 'variety']:
        value = props.get(numeric_field or category_field)
        if value is not None:
            h3_bins[h3_id][category_value]['values'].append(value)
    elif numeric_field:
        value = props.get(numeric_field)
        if value is not None:
            try:
                val = float(value)
                h3_bins[h3_id][category_value][stats].append(val)
            except ValueError:
                pass
