class UPGMATreeBuilder:
    """Unweighted Pair Group Method with Arithmetic Mean (UPGMA) phylogenetic tree builder."""
    def build_tree(self, labels: list[str], distance_matrix: list[list[float]]) -> dict:
        n = len(labels)
        clusters = {i: [i] for i in range(n)}
        node_names = {i: labels[i] for i in range(n)}
        dm = {i: {j: distance_matrix[i][j] for j in range(n)} for i in range(n)}
        merges = []
        next_id = n

        while len(clusters) > 1:
            keys = list(clusters.keys())
            min_d = float('inf')
            best_pair = None

            for i in range(len(keys)):
                for j in range(i + 1, len(keys)):
                    k1, k2 = keys[i], keys[j]
                    if dm[k1][k2] < min_d:
                        min_d = dm[k1][k2]
                        best_pair = (k1, k2)

            k1, k2 = best_pair
            new_cluster = clusters[k1] + clusters[k2]
            new_name = f"({node_names[k1]}:{round(min_d/2.0, 2)},{node_names[k2]}:{round(min_d/2.0, 2)})"
            node_names[next_id] = new_name

            merges.append({
                "cluster1": node_names[k1],
                "cluster2": node_names[k2],
                "distance": round(min_d, 4),
                "branch_height": round(min_d / 2.0, 4)
            })

            # Calculate average distances to new cluster
            new_row = {}
            for other in clusters:
                if other not in (k1, k2):
                    d = sum(distance_matrix[a][b] for a in new_cluster for b in clusters[other]) / (len(new_cluster) * len(clusters[other]))
                    new_row[other] = d
                    dm[other][next_id] = d

            dm[next_id] = new_row
            del clusters[k1]
            del clusters[k2]
            clusters[next_id] = new_cluster
            next_id += 1

        root_id = next_id - 1
        return {
            "num_taxa": n,
            "merges": merges,
            "newick_tree": node_names[root_id] + ";"
        }
