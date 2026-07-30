import sys, statistics
depths = []
with open(sys.argv[1]) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            depths.append(int(parts[2]))
if not depths:
    print("No depth data."); sys.exit(1)
print(f"\n=== Read Depth Summary (chr20) ===")
print(f"Total positions:   {len(depths):,}")
print(f"Mean depth:        {statistics.mean(depths):.2f}x")
print(f"Median depth:      {statistics.median(depths):.2f}x")
print(f"Min / Max:         {min(depths)} / {max(depths)}")
print(f"% positions >=10x: {sum(1 for d in depths if d>=10)/len(depths)*100:.1f}%")
print(f"% positions >=20x: {sum(1 for d in depths if d>=20)/len(depths)*100:.1f}%")
print(f"% positions >=30x: {sum(1 for d in depths if d>=30)/len(depths)*100:.1f}%")
