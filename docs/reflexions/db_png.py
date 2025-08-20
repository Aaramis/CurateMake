from graphviz import Digraph

def generate_updated_schema_png(output_file="its_db_schema.png"):
    """
    Génère un schéma de base de données mise à jour selon la nouvelle documentation.
    """
    # Créer un diagramme
    dot = Digraph(comment="ITS Pipeline Database Schema", format="png")
    
    # --- 1. Configuration globale ---
    dot.attr(rankdir="TB",
             splines="polyline",
             bgcolor="white",
             fontname="Helvetica",
             nodesep="1.5",
             ranksep="2.5",
             pad="0.5")
    
    dot.attr('node', shape="plaintext", fontname="Helvetica")
    dot.attr('edge', fontname="Helvetica", fontsize="13", fontcolor="#2c3e50",
             labeldistance="3.0", labelfontsize="14")

    # --- 2. Palette ---
    colors = {
        "primary": "#2c3e50",
        "bg_table": "#ffffff",
        "pk_bg": "#e8f4fd",
        "fk_bg": "#f0f8ff",
        "edge": "#34495e",
        "border": "#bdc3c7",
        "self_ref": "#7f8c8d"
    }
    
    # --- 3. Tables mises à jour ---
    tables = {
        "Sources": [
            ("source_ID", "PK"),
            ("date_download", ""),
            ("name", "")
        ],
        "Sequences": [
            ("seq_ID", "PK"),
            ("source_ID", "FK"),
            ("original_ID", ""),
            ("sequence", ""),
            ("length", ""),
            ("percentage_N", ""),  # %N renommé pour éviter le %
            ("barcode", ""),
            ("barcode_status", ""), 
            ("parent_seq_ID", "FK") 
        ],
        "Taxonomy": [
            ("tax_ID", "PK"),
            ("kingdom", ""),
            ("phylum", ""),
            ("class", ""),
            ("order", ""),
            ("family", ""),
            ("genus", ""),
            ("species", ""),
            ("parent_tax_ID", "FK")  # Auto-référence hiérarchique
        ],
        "Seq_Tax": [
            ("seq_ID", "FK,PK"),
            ("tax_ID", "FK,PK"),
            ("backbone", ""),  # Ancien source_name -> backbone
            ("status", "")
        ]
    }
    
    # --- 4. Construction des tables ---
    for table_name, fields in tables.items():
        label = f'<<TABLE BORDER="2" CELLBORDER="1" CELLSPACING="0" CELLPADDING="8" COLOR="{colors["border"]}">'
        label += f'<TR><TD BGCOLOR="{colors["primary"]}" ALIGN="CENTER" COLSPAN="2">'
        label += f'<FONT COLOR="white" POINT-SIZE="14"><B>{table_name}</B></FONT></TD></TR>'
        
        for i, (field_name, field_type) in enumerate(fields):
            port_name = f"f{i}"
            
            # Gestion des types de champs
            if "PK" in field_type and "FK" in field_type:
                field_style = f'<I><B>{field_name}</B></I>'  # Italique + Gras
                indicator = '<FONT COLOR="#9b59b6"><B>[FK,PK]</B></FONT>'
                bgcolor = colors["fk_bg"]
            elif "PK" in field_type:
                field_style = f'<B>{field_name}</B>'
                indicator = '<FONT COLOR="#e74c3c"><B>[PK]</B></FONT>'
                bgcolor = colors["pk_bg"]
            elif "FK" in field_type:
                field_style = f'<I>{field_name}</I>'
                indicator = '<FONT COLOR="#3498db"><B>[FK]</B></FONT>'
                bgcolor = colors["fk_bg"]
            else:
                field_style = field_name
                indicator = ""
                bgcolor = colors["bg_table"]
            
            label += (f'<TR>'
                    f'<TD PORT="{port_name}" BGCOLOR="{bgcolor}" ALIGN="LEFT">{field_style}</TD>'
                    f'<TD PORT="{port_name}_pk" BGCOLOR="{bgcolor}" ALIGN="RIGHT">{indicator}</TD>'
                    f'</TR>')
        label += "</TABLE>>"
        dot.node(table_name, label=label)
    
    # --- 5. Relations 1:N classiques ---
    # Sources -> Sequences
    dot.edge("Sources:f0_pk:e", "Sequences:f1",
             taillabel="1", headlabel="N",
             dir="both", arrowhead="crow", arrowtail="tee",
             penwidth="2.5", color=colors["edge"])
    
    # Sequences -> Seq_Tax
    dot.edge("Sequences:f0_pk:e", "Seq_Tax:f0",
             taillabel="1", headlabel="N",
             dir="both", arrowhead="crow", arrowtail="tee",
             penwidth="2.5", color=colors["edge"])
    
    # Taxonomy -> Seq_Tax
    dot.edge("Taxonomy:f0:w", "Seq_Tax:f1_pk:e",
             taillabel="1", headlabel="N",
             dir="both", arrowhead="crow", arrowtail="tee",
             penwidth="2.5", color=colors["edge"])
    
    # --- 6. Relations hiérarchiques (auto-références) ---
    # Sequences.parent_seq_ID -> Sequences.seq_ID
    dot.edge("Sequences:f8_pk:e", "Sequences:f0_pk:e",
             fontsize="11", fontcolor=colors["self_ref"],
             style="dashed", color=colors["self_ref"],
             constraint="false", dir="both", 
             arrowhead="crow", arrowtail="tee")
    
    # Taxonomy.parent_tax_ID -> Taxonomy.tax_ID  
    dot.edge("Taxonomy:f8_pk:e", "Taxonomy:f0_pk:e",
             fontsize="11", fontcolor=colors["self_ref"],
             style="dashed", color=colors["self_ref"],
             constraint="false", dir="both",
             arrowhead="crow", arrowtail="tee")
    
    # --- 7. Titre du diagramme ---
    title = ('<FONT POINT-SIZE="20"><B>Schéma de la base de données - Pipeline ITS (Mise à jour)</B></FONT><BR/>')
    dot.attr(label=f'<{title}>', labelloc="t", fontcolor=colors["primary"])
    
    # --- 8. Légende enrichie ---
    legend = ('<<TABLE BORDER="1" CELLBORDER="0" CELLSPACING="0" CELLPADDING="6" '
              f'COLOR="{colors["border"]}" BGCOLOR="white">'
              f'<TR><TD BGCOLOR="{colors["primary"]}" ALIGN="CENTER">'
              '<FONT COLOR="white" POINT-SIZE="12"><B>Légende</B></FONT></TD></TR>'
              f'<TR><TD BGCOLOR="{colors["pk_bg"]}" ALIGN="LEFT">'
              '<B>Champ</B> <FONT COLOR="#e74c3c"><B>[PK]</B></FONT> = Clé primaire</TD></TR>'
              f'<TR><TD BGCOLOR="{colors["fk_bg"]}" ALIGN="LEFT">'
              '<I>Champ</I> <FONT COLOR="#3498db"><B>[FK]</B></FONT> = Clé étrangère</TD></TR>'
              f'<TR><TD BGCOLOR="{colors["fk_bg"]}" ALIGN="LEFT">'
              '<I><B>Champ</B></I> <FONT COLOR="#9b59b6"><B>[FK,PK]</B></FONT> = Clé composite</TD></TR>'
              '<TR><TD BGCOLOR="white" ALIGN="LEFT">──→ Relation 1:N standard</TD></TR>'
              '<TR><TD BGCOLOR="white" ALIGN="LEFT">┅┅→ Auto-référence hiérarchique</TD></TR>'
              '</TABLE>>')
    dot.node("Legend", label=legend, shape="plaintext")
    
    # --- 9. Organisation spatiale ---
    # Tables principales sur la même ligne
    with dot.subgraph() as s:
        s.attr(rank='same')
        s.node('Sources')
        s.node('Sequences')
        s.node('Seq_Tax')
        s.node('Taxonomy')

    # Arêtes invisibles pour forcer l'ordre horizontal
    dot.edge('Sources', 'Sequences', style='invis')
    dot.edge('Sequences', 'Seq_Tax', style='invis')
    dot.edge('Seq_Tax', 'Taxonomy', style='invis')

    # Légende en bas
    with dot.subgraph() as s:
        s.attr(rank='sink')
        s.node('Legend')
    
    # --- 10. Rendu ---
    dot.render(output_file.replace('.png', ''), view=False, cleanup=True)
    print(f"✅ Schéma mis à jour sauvegardé sous : {output_file}")

if __name__ == "__main__":
    generate_updated_schema_png()