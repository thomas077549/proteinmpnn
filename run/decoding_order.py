"""Serialize sampled decoding orders without changing FASTA sequence content."""


def decoding_fasta_record(record, fasta_chains):
    """List one-based FASTA positions in decoding order under the exact header.

    FASTA chains are alphabetically sorted by the sampler. Number positions
    continuously across those chains, without counting '/' separators. Only
    valid, designed positions are listed; fixed/context residues are omitted.
    """
    positions = {}
    for chain in sorted(fasta_chains):
        for residue in record["residues"]:
            if residue["chain"] == chain:
                positions[residue["model_index"]] = len(positions) + 1
    if len(positions) != len(record["sequence"].replace("/", "")):
        raise ValueError("Decoding-order chain mapping does not match FASTA length")
    order = " ".join(str(positions[index]) for index in record["designed_decoding_order"])
    return f'>{record["fasta_header"]}\n{order}\n'


def decoding_record(protein, chain_letters, order, valid_mask, design_mask,
                    record_index, header, sequence, temperature, noise, model):
    """Indices are model offsets; chain positions are parsed, not PDB author IDs."""
    residues = []
    for chain in chain_letters:
        for position in range(1, len(protein[f"seq_chain_{chain}"]) + 1):
            index = len(residues)
            residues.append({
                "model_index": index, "chain": chain, "chain_position": position,
                "valid": bool(valid_mask[index]), "designed": bool(design_mask[index]),
            })
    if sorted(order) != list(range(len(residues))):
        raise ValueError("Decoding order is not a permutation of the parsed residues")
    return {
        "schema_version": 1,
        "record_index": record_index,
        "fasta_header": header,
        "sequence": sequence,
        "temperature": float(temperature),
        "backbone_noise": float(noise),
        "model_name": model,
        "decoding_order": order,
        "residues": residues,
        "designed_decoding_order": [
            index for index in order
            if residues[index]["valid"] and residues[index]["designed"]
        ],
    }
