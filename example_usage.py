from client import UPGMATreeBuilder

def main():
    print("=== UPGMA Phylogenetic Tree Builder ===")
    builder = UPGMATreeBuilder()
    taxa = ["Human", "Chimp", "Gorilla", "Orangutan"]
    d_mat = [
        [0.0, 1.2, 2.4, 4.8],
        [1.2, 0.0, 2.4, 4.8],
        [2.4, 2.4, 0.0, 4.8],
        [4.8, 4.8, 4.8, 0.0]
    ]

    res = builder.build_tree(taxa, d_mat)
    print("Phylogenetic Tree Newick:", res["newick_tree"])
    assert len(res["merges"]) == 3
    print("UPGMA Tree Builder verified successfully!")

if __name__ == "__main__":
    main()
