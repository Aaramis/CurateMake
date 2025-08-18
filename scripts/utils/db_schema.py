from graphviz import Digraph

# Créer un diagramme entité-relation (ERD)
dot = Digraph(comment="Bioinformatics Pipeline DB Schema", format="png")
dot.attr(rankdir="LR", size="8")

# Définir les tables comme des boîtes
tables = {
    "sources": ["source_id (PK)", "name", "version", "url", "date_download"],
    "sequences": ["seq_id (PK, UUID)", "source_id (FK)", "original_id", "sequence", "length", "md5", "status", "created_at"],
    "taxonomy": ["tax_id (PK)", "kingdom", "phylum", "class", "order", "family", "genus", "species", "scientific_name"],
    "sequence_taxonomy": ["seq_id (FK)", "tax_id (FK)", "method", "confidence", "PRIMARY KEY(seq_id, tax_id, method)"],
    "transformations": ["transform_id (PK)", "seq_id (FK)", "step", "parameters (JSONB)", "tool_version", "created_at"],
    "datasets": ["dataset_id (PK)", "name", "description", "created_at"],
    "dataset_sequences": ["dataset_id (FK)", "seq_id (FK)", "PRIMARY KEY(dataset_id, seq_id)"]
}

# Ajouter les tables au graphe
for table, fields in tables.items():
    label = f"<<TABLE BORDER='1' CELLBORDER='1' CELLSPACING='0'>"
    label += f"<TR><TD BGCOLOR='lightblue'><B>{table}</B></TD></TR>"
    for field in fields:
        label += f"<TR><TD ALIGN='LEFT'>{field}</TD></TR>"
    label += "</TABLE>>"
    dot.node(table, label=label, shape="plaintext")

# Relations
dot.edge("sources", "sequences", label="1:N")
dot.edge("sequences", "sequence_taxonomy", label="1:N")
dot.edge("taxonomy", "sequence_taxonomy", label="1:N")
dot.edge("sequences", "transformations", label="1:N")
dot.edge("datasets", "dataset_sequences", label="1:N")
dot.edge("sequences", "dataset_sequences", label="1:N")

# Rendre l'image
output_path = "/home/agardette/augmentation_taxonomie/5.mvp_snakemake/bioinfo_pipeline_db_schema"
dot.render(output_path, cleanup=True)

output_path + ".png"
