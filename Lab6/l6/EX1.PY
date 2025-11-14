import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def fetch_dna_sequence(min_length=1000, max_length=3000):
    length = random.randint(min_length, max_length)
    sequence = generate_random_sequence(length)
    print(f"Generated sequence of length: {length} bp")

    print("=" * 60)
    print("DNA SEQUENCE")
    print("=" * 60)
    print(sequence)
    print("=" * 60 + "\n")

    return sequence


def generate_random_sequence(length):
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(bases) for _ in range(length))


def generate_fragments(sequence, num_fragments=10, min_size=100, max_size=3000):
    fragments = []
    seq_length = len(sequence)

    for i in range(num_fragments):
        fragment_size = random.randint(min_size, min(max_size, seq_length))

        start_pos = random.randint(0, seq_length - fragment_size)

        fragment = sequence[start_pos:start_pos + fragment_size]
        fragments.append(fragment)

    return fragments


def calculate_migration_distance(fragment_length, max_length=3000, gel_length=10):
    migration = gel_length * (1 - np.log(fragment_length) / np.log(max_length))
    return migration


def visualize_gel_electrophoresis(fragments, title="Gel Electrophoresis Simulation"):
    fig, ax = plt.subplots(figsize=(4, 10), facecolor='black')
    ax.set_facecolor('black')

    gel_width = 3
    gel_height = 10

    gel_box = patches.Rectangle((0.5, 0), gel_width, gel_height,
                                linewidth=0,
                                facecolor='#2a2a2a')
    ax.add_patch(gel_box)

    max_fragment_length = max(len(f) for f in fragments)

    markers = [3000, 1500, 500]
    for marker_size in markers:
        marker_distance = calculate_migration_distance(marker_size, max_fragment_length, gel_length=gel_height - 1)
        marker_y = gel_height - 0.5 - marker_distance
        ax.text(0.2, marker_y, f'{marker_size} bp',
                ha='right', va='center', fontsize=12,
                fontweight='bold', color='white')
        ax.plot([0.5, 0.5 + gel_width], [marker_y, marker_y],
                'r--', alpha=0.2, linewidth=0.5)

    lane_x_start = 0.5 + gel_width * 0.15
    lane_x_width = gel_width * 0.7

    # Draw DNA bands for each fragment in the SAME lane
    for fragment in fragments:
        frag_length = len(fragment)

        migration_distance = calculate_migration_distance(frag_length, max_fragment_length, gel_length=gel_height - 1)

        band_y = gel_height - 0.5 - migration_distance

        band_height = 0.12
        band = patches.Rectangle((lane_x_start, band_y - band_height / 2),
                                 lane_x_width, band_height,
                                 linewidth=0,
                                 facecolor='white',
                                 alpha=0.95)
        ax.add_patch(band)

        ax.text(0.5 + gel_width + 0.1, band_y, f'{frag_length} bp',
                ha='left', va='center', fontsize=10,
                color='white', fontweight='normal')

    ax.set_xlim(0, gel_width + 2)
    ax.set_ylim(-0.5, gel_height + 0.5)
    ax.set_aspect('equal')
    ax.axis('off')

    plt.tight_layout(pad=0.1)
    plt.savefig('gel_electrophoresis.png', dpi=300, bbox_inches='tight', facecolor='black')
    print("\nGel electrophoresis visualization saved as 'gel_electrophoresis.png'")

    return fig


def print_fragment_info(fragments):
    """Print information about the fragments"""
    print("\n" + "=" * 60)
    print("DNA FRAGMENT ANALYSIS")
    print("=" * 60)

    for i, fragment in enumerate(fragments, 1):
        print(f"Fragment {i}: {len(fragment)} bp")

    print("=" * 60)
    sorted_fragments = sorted(enumerate(fragments, 1), key=lambda x: len(x[1]), reverse=True)
    print("\nFragments sorted by size (largest to smallest):")
    for lane, fragment in sorted_fragments:
        print(f"  Fragment {lane}: {len(fragment)} bp")
    print("=" * 60)


def main():
    """Main execution function"""
    print("=" * 60)

    dna_sequence = fetch_dna_sequence(min_length=1000, max_length=3000)

    print("Generating 10 random DNA fragments...")
    fragments = generate_fragments(dna_sequence, num_fragments=10,
                                   min_size=100, max_size=3000)

    print(f"Generated {len(fragments)} fragments")

    print_fragment_info(fragments)

    visualize_gel_electrophoresis(fragments)



if __name__ == "__main__":
    main()