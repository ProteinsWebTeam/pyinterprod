create schema interpro;

create table interpro.protein2go
(
    protein_acc varchar(15) not null,
    term_id varchar(10) not null,
    ref_db_code varchar(10) not null,
    ref_db_id varchar(60) not null
);

create index protein2go_protein_idx
    on protein2go (protein_acc);

create index protein2go_term_idx
    on protein2go (term_id);

create table interpro.publication
(
    id varchar(25) not null
        constraint publication_pkey
            primary key,
    title varchar(1500) not null,
    published date not null
);

create table interpro.term
(
    id varchar(10) not null
        constraint term_pkey
            primary key,
    name varchar(200) not null,
    category varchar(25) not null,
    num_constraints integer not null,
    is_obsolete boolean not null,
    definition varchar not null,
    replaced_id varchar(10)
);

create table interpro.lineage
(
    child_id integer not null,
    parent_id integer not null,
    parent_rank varchar(255) not null
);

create unique index lineage_child_parent_uidx
    on lineage (child_id, parent_id);

create index lineage_child_rank_idx
    on lineage (child_id, parent_rank);

create table interpro.protein
(
    accession varchar(15) not null
        constraint protein_pkey
            primary key,
    identifier varchar(16) not null,
    length integer not null,
    taxon_id integer not null,
    is_fragment boolean not null,
    is_reviewed boolean not null
);

create unique index protein_identifier_uidx
    on protein (identifier);

create index protein_reviewed_idx
    on protein (is_reviewed);

create table interpro.match
(
    protein_acc varchar(15) not null,
    signature_acc varchar(25) not null,
    database_id integer not null,
    fragments text not null
);

create index match_protein_idx
    on match (protein_acc);

create table interpro.protein2name
(
    protein_acc varchar(15) not null
        constraint protein2name_pkey
            primary key,
    name_id integer not null
);

create table interpro.protein_name
(
    name_id integer not null
        constraint protein_name_pkey
            primary key,
    text text not null
);

create table interpro.signature2protein
(
    signature_acc varchar(25) not null,
    protein_acc varchar(15) not null,
    is_reviewed boolean not null,
    taxon_left_num integer not null,
    name_id integer not null
);

create index signature2protein_protein_idx
    on signature2protein (protein_acc);

create index signature2protein_signature_idx
    on signature2protein (signature_acc);

create table interpro.taxon
(
    id integer not null
        constraint taxon_id_pkey
            primary key,
    name varchar(255) not null,
    rank varchar(50) not null,
    left_number integer not null,
    right_number integer not null,
    parent_id integer,
    lineage text not null
);

create index taxon_left_number_idx
    on taxon (left_number);

create table interpro.database
(
    id serial not null
        constraint database_pkey
            primary key,
    name varchar(50) not null,
    name_long varchar(50) not null,
    version varchar(20) not null,
    updated date not null
);

create unique index database_name_idx
    on database (name);

create table interpro.signature
(
    accession varchar(25) not null
        constraint signature_pkey
            primary key,
    database_id integer not null,
    name varchar(100) not null,
    description varchar(400),
    type varchar(25) not null,
    abstract text,
    num_sequences integer not null,
    num_complete_sequences integer not null,
    num_residues bigint not null
);

create index signature_database_idx
    on signature (database_id);

create table interpro.protein_similarity
(
    comment_id integer not null,
    comment_text text not null,
    protein_acc varchar(15) not null
);

create index protein_similarity_protein_idx
    on protein_similarity (protein_acc);

create index protein_similarity_comment_idx
    on protein_similarity (comment_id);

create table interpro.comparison
(
    signature_acc_1 varchar(25) not null,
    signature_acc_2 varchar(25) not null,
    collocations integer not null,
    "overlaps" integer not null
);

create index comparison_idx
    on comparison (signature_acc_1, signature_acc_2);

create table interpro.prediction
(
    signature_acc_1 varchar(25) not null,
    signature_acc_2 varchar(25) not null,
    collocations integer not null,
    protein_overlaps integer not null,
    residue_overlaps integer not null
);

create index prediction_idx
    on prediction (signature_acc_1);
