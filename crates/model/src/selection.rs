//! A small node selection language for documents.
//!
//! A selection is a boolean expression over three node kinds: parts, atoms,
//! and bonds. Each expression evaluates to one set per kind. `and`, `or`, and
//! `not` apply the set algebra of each kind separately. A predicate that names
//! atoms leaves the bond and part sets empty, and the reverse holds too.

use std::collections::BTreeSet;

use crate::document::Document;
use crate::element::Element;
use crate::error::ModelError;
use crate::part::Part;

/// A parsed node selection expression.
///
/// Build one with [`parse_selection`] and run it against a document with
/// [`Selection::evaluate`].
#[derive(Clone, Debug, PartialEq)]
pub struct Selection {
    expr: Expr,
}

/// The node sets a selection selects.
///
/// `parts` holds part indices. `atoms` and `bonds` hold `(part_index,
/// node_index)` pairs. Every list is ascending and free of duplicates.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SelectionResult {
    /// The selected part indices.
    pub parts: Vec<usize>,
    /// The selected `(part_index, atom_index)` pairs.
    pub atoms: Vec<(usize, usize)>,
    /// The selected `(part_index, bond_index)` pairs.
    pub bonds: Vec<(usize, usize)>,
}

/// Parses a selection expression.
///
/// Keywords are case-insensitive and whitespace does not matter. The function
/// returns [`ModelError::SelectionParse`] when the text is not valid.
pub fn parse_selection(text: &str) -> Result<Selection, ModelError> {
    let tokens = lex(text)?;
    let mut parser = Parser {
        tokens,
        index: 0,
        end: text.len(),
    };
    let expr = parser.parse_or()?;
    let (token, position) = parser
        .peek()
        .ok_or_else(|| selection_error(text.len(), "unexpected end of selection"))?;
    if !matches!(token, Token::Eof) {
        return Err(selection_error(*position, "unexpected trailing input"));
    }
    Ok(Selection { expr })
}

impl Selection {
    /// Evaluates the selection against a document.
    ///
    /// The result is deterministic. It lists every selected node in ascending
    /// order without duplicates.
    pub fn evaluate(&self, document: &Document) -> SelectionResult {
        let full = full_sets(document);
        let sets = evaluate_expr(&self.expr, document, &full);
        SelectionResult {
            parts: sets.parts.into_iter().collect(),
            atoms: sets.atoms.into_iter().collect(),
            bonds: sets.bonds.into_iter().collect(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Sets {
    parts: BTreeSet<usize>,
    atoms: BTreeSet<(usize, usize)>,
    bonds: BTreeSet<(usize, usize)>,
}

#[derive(Clone, Debug, PartialEq)]
enum Expr {
    All,
    None,
    Kind(NodeKind),
    Element(Element),
    Charge(Cmp, f64),
    Degree(Cmp, f64),
    AtomIndex(Cmp, f64),
    Order(Cmp, f64),
    BondLength(Cmp, f64),
    BondIndex(Cmp, f64),
    PartIndex(Cmp, f64),
    Name(String),
    Within { distance: f64, atom: usize },
    Not(Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum NodeKind {
    Atom,
    Bond,
    Part,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
}

impl Cmp {
    fn test(self, left: f64, right: f64) -> bool {
        match self {
            Self::Lt => left < right,
            Self::Le => left <= right,
            Self::Gt => left > right,
            Self::Ge => left >= right,
            Self::Eq => left == right,
            Self::Ne => left != right,
        }
    }
}

fn selection_error(position: usize, message: impl Into<String>) -> ModelError {
    ModelError::SelectionParse {
        position,
        message: message.into(),
    }
}

fn full_sets(document: &Document) -> Sets {
    let mut sets = Sets::default();
    for (part_index, part) in document.parts.iter().enumerate() {
        sets.parts.insert(part_index);
        for atom_index in 0..part.atom_count() {
            sets.atoms.insert((part_index, atom_index));
        }
        for bond_index in 0..part.bond_count() {
            sets.bonds.insert((part_index, bond_index));
        }
    }
    sets
}

fn evaluate_expr(expr: &Expr, document: &Document, full: &Sets) -> Sets {
    match expr {
        Expr::All => full.clone(),
        Expr::None => Sets::default(),
        Expr::Kind(NodeKind::Atom) => Sets {
            atoms: full.atoms.clone(),
            ..Sets::default()
        },
        Expr::Kind(NodeKind::Bond) => Sets {
            bonds: full.bonds.clone(),
            ..Sets::default()
        },
        Expr::Kind(NodeKind::Part) => Sets {
            parts: full.parts.clone(),
            ..Sets::default()
        },
        Expr::Element(element) => atoms_where(document, full, |index, part| {
            part.topology.element(index) == Some(*element)
        }),
        Expr::Charge(cmp, value) => atoms_where(document, full, |index, part| {
            part.topology
                .charge_c(index)
                .is_some_and(|charge| cmp.test(charge, *value))
        }),
        Expr::Degree(cmp, value) => atoms_where(document, full, |index, part| {
            cmp.test(degree(part, index) as f64, *value)
        }),
        Expr::AtomIndex(cmp, value) => {
            atoms_where(document, full, |index, _| cmp.test(index as f64, *value))
        }
        Expr::Order(cmp, value) => bonds_where(document, full, |index, part| {
            part.topology
                .bond_order(index)
                .is_some_and(|order| cmp.test(order as f64, *value))
        }),
        Expr::BondLength(cmp, value) => bonds_where(document, full, |index, part| {
            bond_length(part, index).is_some_and(|length| cmp.test(length, *value))
        }),
        Expr::BondIndex(cmp, value) => {
            bonds_where(document, full, |index, _| cmp.test(index as f64, *value))
        }
        Expr::PartIndex(cmp, value) => {
            parts_where(document, full, |index, _| cmp.test(index as f64, *value))
        }
        Expr::Name(name) => parts_where(document, full, |_, part| part.name == *name),
        Expr::Within { distance, atom } => within(document, *distance, *atom),
        Expr::Not(inner) => complement(&evaluate_expr(inner, document, full), full),
        Expr::And(left, right) => intersect(
            &evaluate_expr(left, document, full),
            &evaluate_expr(right, document, full),
        ),
        Expr::Or(left, right) => union(
            &evaluate_expr(left, document, full),
            &evaluate_expr(right, document, full),
        ),
    }
}

fn atoms_where(document: &Document, full: &Sets, predicate: impl Fn(usize, &Part) -> bool) -> Sets {
    let mut sets = Sets::default();
    for &(part_index, atom_index) in &full.atoms {
        if let Some(part) = document.part(part_index) {
            if predicate(atom_index, part) {
                sets.atoms.insert((part_index, atom_index));
            }
        }
    }
    sets
}

fn bonds_where(document: &Document, full: &Sets, predicate: impl Fn(usize, &Part) -> bool) -> Sets {
    let mut sets = Sets::default();
    for &(part_index, bond_index) in &full.bonds {
        if let Some(part) = document.part(part_index) {
            if predicate(bond_index, part) {
                sets.bonds.insert((part_index, bond_index));
            }
        }
    }
    sets
}

fn parts_where(document: &Document, full: &Sets, predicate: impl Fn(usize, &Part) -> bool) -> Sets {
    let mut sets = Sets::default();
    for &part_index in &full.parts {
        if let Some(part) = document.part(part_index) {
            if predicate(part_index, part) {
                sets.parts.insert(part_index);
            }
        }
    }
    sets
}

fn degree(part: &Part, atom_index: usize) -> usize {
    (0..part.bond_count())
        .filter(|&bond_index| {
            part.topology
                .bond(bond_index)
                .is_some_and(|bond| bond.u as usize == atom_index || bond.v as usize == atom_index)
        })
        .count()
}

fn bond_length(part: &Part, bond_index: usize) -> Option<f64> {
    let bond = part.topology.bond(bond_index)?;
    let u = part.topology.position_m(bond.u as usize)?;
    let v = part.topology.position_m(bond.v as usize)?;
    Some(distance(u, v))
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn within(document: &Document, max_distance: f64, atom_index: usize) -> Sets {
    let mut sets = Sets::default();
    for (part_index, part) in document.parts.iter().enumerate() {
        let Some(centre) = part.topology.position_m(atom_index) else {
            continue;
        };
        for candidate in 0..part.atom_count() {
            let Some(position) = part.topology.position_m(candidate) else {
                continue;
            };
            if distance(centre, position) <= max_distance {
                sets.atoms.insert((part_index, candidate));
            }
        }
    }
    sets
}

fn complement(inner: &Sets, full: &Sets) -> Sets {
    Sets {
        parts: full.parts.difference(&inner.parts).copied().collect(),
        atoms: full.atoms.difference(&inner.atoms).copied().collect(),
        bonds: full.bonds.difference(&inner.bonds).copied().collect(),
    }
}

fn intersect(left: &Sets, right: &Sets) -> Sets {
    Sets {
        parts: left.parts.intersection(&right.parts).copied().collect(),
        atoms: left.atoms.intersection(&right.atoms).copied().collect(),
        bonds: left.bonds.intersection(&right.bonds).copied().collect(),
    }
}

fn union(left: &Sets, right: &Sets) -> Sets {
    Sets {
        parts: left.parts.union(&right.parts).copied().collect(),
        atoms: left.atoms.union(&right.atoms).copied().collect(),
        bonds: left.bonds.union(&right.bonds).copied().collect(),
    }
}

#[derive(Clone, Debug, PartialEq)]
enum Token {
    Ident(String),
    Number(f64),
    Str(String),
    Cmp(Cmp),
    LParen,
    RParen,
    LBracket,
    RBracket,
    Eof,
}

fn lex(text: &str) -> Result<Vec<(Token, usize)>, ModelError> {
    let mut lexer = Lexer::new(text);
    let mut tokens = Vec::new();
    while let Some((token, position)) = lexer.next_token()? {
        tokens.push((token, position));
    }
    tokens.push((Token::Eof, text.len()));
    Ok(tokens)
}

struct Lexer<'a> {
    text: &'a str,
    pos: usize,
}

impl<'a> Lexer<'a> {
    fn new(text: &'a str) -> Self {
        Self { text, pos: 0 }
    }

    fn peek_char(&self) -> Option<char> {
        self.text[self.pos..].chars().next()
    }

    fn bump(&mut self) -> Option<char> {
        let ch = self.peek_char()?;
        self.pos += ch.len_utf8();
        Some(ch)
    }

    fn skip_whitespace(&mut self) {
        while self.peek_char().is_some_and(char::is_whitespace) {
            self.bump();
        }
    }

    fn next_token(&mut self) -> Result<Option<(Token, usize)>, ModelError> {
        self.skip_whitespace();
        let start = self.pos;
        let Some(ch) = self.peek_char() else {
            return Ok(None);
        };
        let token = match ch {
            '(' => {
                self.bump();
                Token::LParen
            }
            ')' => {
                self.bump();
                Token::RParen
            }
            '[' => {
                self.bump();
                Token::LBracket
            }
            ']' => {
                self.bump();
                Token::RBracket
            }
            '"' => Token::Str(self.lex_string(start)?),
            '<' => {
                self.bump();
                if self.peek_char() == Some('=') {
                    self.bump();
                    Token::Cmp(Cmp::Le)
                } else {
                    Token::Cmp(Cmp::Lt)
                }
            }
            '>' => {
                self.bump();
                if self.peek_char() == Some('=') {
                    self.bump();
                    Token::Cmp(Cmp::Ge)
                } else {
                    Token::Cmp(Cmp::Gt)
                }
            }
            '=' => {
                self.bump();
                if self.peek_char() == Some('=') {
                    self.bump();
                    Token::Cmp(Cmp::Eq)
                } else {
                    return Err(selection_error(start, "expected `==`"));
                }
            }
            '!' => {
                self.bump();
                if self.peek_char() == Some('=') {
                    self.bump();
                    Token::Cmp(Cmp::Ne)
                } else {
                    return Err(selection_error(start, "expected `!=`"));
                }
            }
            _ if ch.is_ascii_digit() || ch == '.' => Token::Number(self.lex_number(start)?),
            _ if ch.is_alphabetic() || ch == '_' => Token::Ident(self.lex_ident()),
            _ => return Err(selection_error(start, "unexpected character")),
        };
        Ok(Some((token, start)))
    }

    fn lex_ident(&mut self) -> String {
        let start = self.pos;
        while self
            .peek_char()
            .is_some_and(|ch| ch.is_alphanumeric() || ch == '_')
        {
            self.bump();
        }
        self.text[start..self.pos].to_string()
    }

    fn lex_number(&mut self, start: usize) -> Result<f64, ModelError> {
        while self
            .peek_char()
            .is_some_and(|ch| ch.is_ascii_digit() || ch == '.')
        {
            self.bump();
        }
        if self.peek_char().is_some_and(|ch| ch == 'e' || ch == 'E') {
            self.bump();
            if self.peek_char().is_some_and(|ch| ch == '+' || ch == '-') {
                self.bump();
            }
            while self.peek_char().is_some_and(|ch| ch.is_ascii_digit()) {
                self.bump();
            }
        }
        let slice = &self.text[start..self.pos];
        slice
            .parse::<f64>()
            .map_err(|_| selection_error(start, "invalid number literal"))
    }

    fn lex_string(&mut self, start: usize) -> Result<String, ModelError> {
        self.bump();
        let mut value = String::new();
        loop {
            let Some(ch) = self.bump() else {
                return Err(selection_error(start, "unterminated string"));
            };
            match ch {
                '"' => return Ok(value),
                '\\' => {
                    let Some(escaped) = self.bump() else {
                        return Err(selection_error(start, "unterminated string"));
                    };
                    match escaped {
                        '"' => value.push('"'),
                        '\\' => value.push('\\'),
                        _ => {
                            return Err(selection_error(
                                self.pos - escaped.len_utf8(),
                                "unknown string escape",
                            ))
                        }
                    }
                }
                _ => value.push(ch),
            }
        }
    }
}

struct Parser {
    tokens: Vec<(Token, usize)>,
    index: usize,
    end: usize,
}

impl Parser {
    fn peek(&self) -> Option<&(Token, usize)> {
        self.tokens.get(self.index)
    }

    fn advance(&mut self) -> Option<(Token, usize)> {
        let token = self.tokens.get(self.index).cloned();
        if token.is_some() {
            self.index += 1;
        }
        token
    }

    fn position(&self) -> usize {
        self.peek().map_or(self.end, |(_, position)| *position)
    }

    fn take_cmp(&mut self, context: &str) -> Result<Cmp, ModelError> {
        match self.advance() {
            Some((Token::Cmp(cmp), _)) => Ok(cmp),
            _ => Err(selection_error(
                self.position(),
                format!("expected a comparison operator after `{context}`"),
            )),
        }
    }

    fn take_number(&mut self) -> Result<f64, ModelError> {
        match self.advance() {
            Some((Token::Number(value), _)) => Ok(value),
            _ => Err(selection_error(self.position(), "expected a number")),
        }
    }

    fn take_string(&mut self) -> Result<String, ModelError> {
        match self.advance() {
            Some((Token::Str(value), _)) => Ok(value),
            _ => Err(selection_error(self.position(), "expected a string")),
        }
    }

    fn take_keyword(&mut self, keyword: &str) -> Result<(), ModelError> {
        if self.peek_keyword(keyword) {
            self.advance();
            Ok(())
        } else {
            Err(selection_error(
                self.position(),
                format!("expected `{keyword}`"),
            ))
        }
    }

    fn peek_keyword(&self, keyword: &str) -> bool {
        matches!(
            self.peek(),
            Some((Token::Ident(name), _)) if name.eq_ignore_ascii_case(keyword)
        )
    }

    fn parse_or(&mut self) -> Result<Expr, ModelError> {
        let mut left = self.parse_and()?;
        while self.peek_keyword("or") {
            self.advance();
            let right = self.parse_and()?;
            left = Expr::Or(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn parse_and(&mut self) -> Result<Expr, ModelError> {
        let mut left = self.parse_factor()?;
        while self.peek_keyword("and") {
            self.advance();
            let right = self.parse_factor()?;
            left = Expr::And(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn parse_factor(&mut self) -> Result<Expr, ModelError> {
        if self.peek_keyword("not") {
            self.advance();
            return Ok(Expr::Not(Box::new(self.parse_factor()?)));
        }
        if matches!(self.peek(), Some((Token::LParen, _))) {
            self.advance();
            let expr = self.parse_or()?;
            match self.advance() {
                Some((Token::RParen, _)) => return Ok(expr),
                _ => return Err(selection_error(self.position(), "expected `)`")),
            }
        }
        self.parse_predicate()
    }

    fn parse_predicate(&mut self) -> Result<Expr, ModelError> {
        let Some((token, position)) = self.advance() else {
            return Err(selection_error(self.end, "expected a predicate"));
        };
        let Token::Ident(name) = token else {
            return Err(selection_error(position, "expected a predicate"));
        };
        let keyword = name.to_ascii_lowercase();
        let expr = match keyword.as_str() {
            "all" => Expr::All,
            "none" => Expr::None,
            "atom" => Expr::Kind(NodeKind::Atom),
            "bond" => Expr::Kind(NodeKind::Bond),
            "part" => Expr::Kind(NodeKind::Part),
            "carbon" => Expr::Element(Element::CARBON),
            "hydrogen" => Expr::Element(Element::HYDROGEN),
            "nitrogen" => Expr::Element(Element::NITROGEN),
            "oxygen" => Expr::Element(Element::OXYGEN),
            "charge" => Expr::Charge(self.take_cmp("charge")?, self.take_number()?),
            "degree" => Expr::Degree(self.take_cmp("degree")?, self.take_number()?),
            "atom_index" => Expr::AtomIndex(self.take_cmp("atom_index")?, self.take_number()?),
            "order" => Expr::Order(self.take_cmp("order")?, self.take_number()?),
            "bond_length" => Expr::BondLength(self.take_cmp("bond_length")?, self.take_number()?),
            "bond_index" => Expr::BondIndex(self.take_cmp("bond_index")?, self.take_number()?),
            "part_index" => Expr::PartIndex(self.take_cmp("part_index")?, self.take_number()?),
            "name" => {
                let cmp = self.take_cmp("name")?;
                if cmp != Cmp::Eq {
                    return Err(selection_error(position, "`name` supports only `==`"));
                }
                Expr::Name(self.take_string()?)
            }
            "within" => {
                let max_distance = self.take_number()?;
                self.take_keyword("of")?;
                self.take_keyword("atom")?;
                match self.advance() {
                    Some((Token::LBracket, _)) => {}
                    _ => return Err(selection_error(self.position(), "expected `[`")),
                }
                let reference_position = self.position();
                let reference = self.take_number()?;
                match self.advance() {
                    Some((Token::RBracket, _)) => {}
                    _ => return Err(selection_error(self.position(), "expected `]`")),
                }
                let atom = atom_reference(reference, reference_position)?;
                Expr::Within {
                    distance: max_distance,
                    atom,
                }
            }
            _ => {
                return Err(selection_error(
                    position,
                    format!("unknown predicate `{name}`"),
                ))
            }
        };
        Ok(expr)
    }
}

fn atom_reference(value: f64, position: usize) -> Result<usize, ModelError> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 {
        return Err(selection_error(
            position,
            "atom reference must be a non-negative integer",
        ));
    }
    Ok(value as usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate as nanocad_model;
    use nanocad_model::{Atom, Bond, BondType, Document, Element, Part, Topology};

    fn chain() -> Part {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(Element::CARBON, [0.0, 0.0, 0.0], 0.3, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [1.0e-10, 0.0, 0.0], 0.0, "C3"));
        topology.add_atom(Atom::new(Element::CARBON, [2.0e-10, 0.0, 0.0], 0.0, "C3"));
        topology
            .add_bond(Bond::new(0, 1, 1, BondType::Single))
            .expect("valid bond");
        topology
            .add_bond(Bond::new(1, 2, 1, BondType::Single))
            .expect("valid bond");
        Part::new("sun", topology)
    }

    fn hydrogen() -> Part {
        let mut topology = Topology::new();
        topology.add_atom(Atom::new(
            Element::HYDROGEN,
            [9.0e-10, 9.0e-10, 9.0e-10],
            0.0,
            "H1",
        ));
        Part::new("gear", topology)
    }

    fn document() -> Document {
        let mut document = Document::new("fixture");
        document.add_part(chain());
        document.add_part(hydrogen());
        document
    }

    fn select(text: &str) -> SelectionResult {
        parse_selection(text)
            .expect("valid selection")
            .evaluate(&document())
    }

    #[test]
    fn carbon_selects_the_three_carbons_and_nothing_else() {
        let result = select("carbon");
        assert_eq!(result.atoms, vec![(0, 0), (0, 1), (0, 2)]);
        assert!(result.parts.is_empty());
        assert!(result.bonds.is_empty());
    }

    #[test]
    fn carbon_and_degree_selects_the_middle_carbon() {
        let result = select("carbon and degree == 2");
        assert_eq!(result.atoms, vec![(0, 1)]);
    }

    #[test]
    fn not_carbon_selects_the_hydrogen() {
        let result = select("not carbon");
        assert_eq!(result.atoms, vec![(1, 0)]);
        assert!(result.atoms.iter().all(|&(part, _)| part == 1));
    }

    #[test]
    fn atom_index_and_charge_filter_atoms() {
        let by_index = select("atom_index < 1");
        assert_eq!(by_index.atoms, vec![(0, 0), (1, 0)]);
        let by_charge = select("charge > 0.0");
        assert_eq!(by_charge.atoms, vec![(0, 0)]);
    }

    #[test]
    fn name_selects_the_named_part() {
        let result = select("name == \"gear\"");
        assert_eq!(result.parts, vec![1]);
        assert!(result.atoms.is_empty());
    }

    #[test]
    fn bond_length_selects_by_length() {
        let result = select("bond_length < 2.0e-10");
        assert_eq!(result.bonds, vec![(0, 0), (0, 1)]);
        assert!(result.atoms.is_empty());
    }

    #[test]
    fn parse_errors_carry_a_position_and_do_not_panic() {
        let error = parse_selection("carbon and").expect_err("must fail");
        match error {
            ModelError::SelectionParse { position, .. } => {
                assert_eq!(position, "carbon and".len());
            }
            other => panic!("unexpected error: {other:?}"),
        }
        let unknown = parse_selection("carbon and frobnicate").expect_err("must fail");
        match unknown {
            ModelError::SelectionParse { position, message } => {
                assert_eq!(position, "carbon and ".len());
                assert!(message.contains("frobnicate"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn all_selects_every_node_and_none_selects_nothing() {
        let result = select("all");
        assert_eq!(result.parts, vec![0, 1]);
        assert_eq!(result.atoms, vec![(0, 0), (0, 1), (0, 2), (1, 0)]);
        assert_eq!(result.bonds, vec![(0, 0), (0, 1)]);
        let empty = select("none");
        assert!(empty.parts.is_empty());
        assert!(empty.atoms.is_empty());
        assert!(empty.bonds.is_empty());
    }

    #[test]
    fn within_and_of_atom_stays_in_one_part() {
        let result = select("within 1.5e-10 of atom [1]");
        assert_eq!(result.atoms, vec![(0, 0), (0, 1), (0, 2)]);
        assert!(result.atoms.iter().all(|&(part, _)| part == 0));
    }

    #[test]
    fn keywords_are_case_insensitive() {
        let result = select("CARBON AND Degree == 2");
        assert_eq!(result.atoms, vec![(0, 1)]);
    }

    #[test]
    fn not_is_the_complement_of_the_document() {
        let result = select("not atom");
        assert!(result.atoms.is_empty());
        assert_eq!(result.parts, vec![0, 1]);
        assert_eq!(result.bonds, vec![(0, 0), (0, 1)]);
    }
}
