use std::collections::BTreeMap;

use crate::error::ModelError;

pub(crate) struct Writer {
    buf: Vec<u8>,
}

impl Writer {
    pub(crate) fn new() -> Self {
        Self { buf: Vec::new() }
    }

    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.buf
    }

    pub(crate) fn write_u8(&mut self, value: u8) {
        self.buf.push(value);
    }

    pub(crate) fn write_u32(&mut self, value: u32) {
        self.buf.extend_from_slice(&value.to_le_bytes());
    }

    pub(crate) fn write_f64(&mut self, value: f64) {
        self.buf.extend_from_slice(&value.to_le_bytes());
    }

    pub(crate) fn write_bytes(&mut self, bytes: &[u8]) {
        self.buf.extend_from_slice(bytes);
    }

    pub(crate) fn write_str(&mut self, value: &str) -> Result<(), ModelError> {
        let len = u32::try_from(value.len()).map_err(|_| ModelError::IntegerOverflow)?;
        self.write_u32(len);
        self.write_bytes(value.as_bytes());
        Ok(())
    }

    pub(crate) fn write_count(&mut self, count: usize) -> Result<(), ModelError> {
        let count = u32::try_from(count).map_err(|_| ModelError::IntegerOverflow)?;
        self.write_u32(count);
        Ok(())
    }

    pub(crate) fn write_map(&mut self, map: &BTreeMap<String, String>) -> Result<(), ModelError> {
        self.write_count(map.len())?;
        for (key, value) in map {
            self.write_str(key)?;
            self.write_str(value)?;
        }
        Ok(())
    }
}

pub(crate) struct Reader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    pub(crate) fn new(buf: &'a [u8]) -> Self {
        Self { buf, pos: 0 }
    }

    pub(crate) fn remaining(&self) -> usize {
        self.buf.len().saturating_sub(self.pos)
    }

    fn read_exact(&mut self, len: usize) -> Result<&'a [u8], ModelError> {
        let end = self
            .pos
            .checked_add(len)
            .ok_or(ModelError::UnexpectedEof { offset: self.pos })?;
        let slice = self
            .buf
            .get(self.pos..end)
            .ok_or(ModelError::UnexpectedEof { offset: self.pos })?;
        self.pos = end;
        Ok(slice)
    }

    pub(crate) fn read_u8(&mut self) -> Result<u8, ModelError> {
        let bytes = self.read_exact(1)?;
        Ok(bytes[0])
    }

    pub(crate) fn read_array<const N: usize>(&mut self) -> Result<[u8; N], ModelError> {
        let bytes = self.read_exact(N)?;
        let mut array = [0u8; N];
        array.copy_from_slice(bytes);
        Ok(array)
    }

    pub(crate) fn read_u32(&mut self) -> Result<u32, ModelError> {
        let bytes = self.read_exact(4)?;
        Ok(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    pub(crate) fn read_f64(&mut self) -> Result<f64, ModelError> {
        let bytes = self.read_exact(8)?;
        Ok(f64::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
        ]))
    }

    pub(crate) fn read_str(&mut self) -> Result<String, ModelError> {
        let len = self.read_u32()? as usize;
        let bytes = self.read_exact(len)?;
        std::str::from_utf8(bytes)
            .map(str::to_owned)
            .map_err(|_| ModelError::InvalidUtf8)
    }

    pub(crate) fn read_count(&mut self) -> Result<usize, ModelError> {
        Ok(self.read_u32()? as usize)
    }

    pub(crate) fn read_map(&mut self) -> Result<BTreeMap<String, String>, ModelError> {
        let count = self.read_count()?;
        let mut map = BTreeMap::new();
        for _ in 0..count {
            let key = self.read_str()?;
            let value = self.read_str()?;
            map.insert(key, value);
        }
        Ok(map)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integers_and_floats_round_trip_through_the_writer() {
        let mut writer = Writer::new();
        writer.write_u8(7);
        writer.write_u32(0xdead_beef);
        writer.write_f64(-1.25e-10);
        let bytes = writer.into_bytes();

        let mut reader = Reader::new(&bytes);
        assert_eq!(reader.read_u8().expect("u8"), 7);
        assert_eq!(reader.read_u32().expect("u32"), 0xdead_beef);
        assert_eq!(reader.read_f64().expect("f64"), -1.25e-10);
        assert_eq!(reader.remaining(), 0);
    }

    #[test]
    fn a_reading_past_the_end_is_an_error_not_a_panic() {
        let bytes = [1u8, 2, 3];
        let mut reader = Reader::new(&bytes);
        assert_eq!(
            reader.read_u32(),
            Err(ModelError::UnexpectedEof { offset: 0 })
        );
    }

    #[test]
    fn a_string_round_trips() {
        let mut writer = Writer::new();
        writer.write_str("diamond").expect("fits");
        let bytes = writer.into_bytes();
        let mut reader = Reader::new(&bytes);
        assert_eq!(reader.read_str().expect("string"), "diamond");
    }

    #[test]
    fn a_map_round_trips_in_key_order() {
        let mut map = BTreeMap::new();
        map.insert("b".to_owned(), "2".to_owned());
        map.insert("a".to_owned(), "1".to_owned());
        let mut writer = Writer::new();
        writer.write_map(&map).expect("fits");
        let bytes = writer.into_bytes();
        let mut reader = Reader::new(&bytes);
        assert_eq!(reader.read_map().expect("map"), map);
    }

    #[test]
    fn an_invalid_utf8_string_is_rejected() {
        let mut writer = Writer::new();
        writer.write_u32(2);
        writer.write_bytes(&[0xff, 0xfe]);
        let bytes = writer.into_bytes();
        let mut reader = Reader::new(&bytes);
        assert_eq!(reader.read_str(), Err(ModelError::InvalidUtf8));
    }
}
