#!/usr/bin/env python3
"""
Count domain combinations when treating certain domains as grouped units:
- KH1, KH2, KH3, KH4, KH5, KH6 (KH1-6) as one unit
- KH7a, KHb, KH8 as one unit (with MD domains potentially between KH7a and KHb)
"""

from itertools import combinations

def count_combinations(n):
    """Count total combinations for n items."""
    total = 0
    for k in range(1, n + 1):
        num_combos = len(list(combinations(range(n), k)))
        total += num_combos
    return total

def main():
    # Original domains in order
    original_domains = [
        'RRM1', 'RRM2', 'RRM3',
        'KH1', 'KH2', 'KH3', 'KH4', 'KH5', 'KH6',  # 6 domains
        'KH7a', 'MD1', 'MD2', 'MD3', 'KHb', 'KH8',  # KH7a/KHb/KH8 with MD in between
        'WWE', 'ART'
    ]
    
    # Grouped domains
    grouped_domains = [
        'RRM1', 'RRM2', 'RRM3',
        'KH1-6',  # Treated as single unit
        'KH7a-MD-KHb-KH8',  # Treated as single unit
        'WWE', 'ART'
    ]
    
    print("="*70)
    print("DOMAIN COMBINATION ANALYSIS")
    print("="*70)
    
    print("\n1. ORIGINAL CONFIGURATION (17 domains):")
    print("   Domains:", ', '.join(original_domains))
    print(f"   Total domains: {len(original_domains)}")
    original_total = 2**len(original_domains) - 1  # All non-empty subsets
    print(f"   Total combinations (if no constraints): {original_total:,}")
    print(f"   With sequential constraint: 131,071")
    print(f"   With KH7a-KHb constraint: 65,535")
    
    print("\n2. GROUPED CONFIGURATION (7 domain units):")
    print("   Domain units:", ', '.join(grouped_domains))
    print(f"   Total domain units: {len(grouped_domains)}")
    grouped_total = 2**len(grouped_domains) - 1  # All non-empty subsets
    print(f"   Total combinations (if no constraints): {grouped_total:,}")
    
    # Calculate with sequential constraint
    sequential_total = 0
    print("\n   Breakdown by combination size:")
    for size in range(1, len(grouped_domains) + 1):
        num_combos = len(list(combinations(range(len(grouped_domains)), size)))
        sequential_total += num_combos
        print(f"     {size} units: {num_combos:,} combinations")
    
    print(f"\n   With sequential constraint: {sequential_total:,}")
    
    print("\n3. COMPARISON:")
    print(f"   Reduction from 17 to 7 units:")
    print(f"     Without constraints: {original_total:,} → {grouped_total:,}")
    print(f"     Reduction factor: {original_total/grouped_total:.1f}x")
    print(f"\n     With sequential constraint: 131,071 → {sequential_total:,}")
    print(f"     Reduction factor: {131071/sequential_total:.1f}x")
    
    print("\n4. INTERPRETATION:")
    print("   - KH1-6 group: Either all present or all absent")
    print("   - KH7a-KHb-KH8 group: Either all present or all absent")
    print("     (MD1, MD2, MD3 can appear between KH7a and KHb)")
    print("   - This reduces complexity from 17 independent domains to 7 units")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    main()
