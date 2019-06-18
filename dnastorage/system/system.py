

SystemFileFormats = {
    # KEY      KEY    LEN  Abbrev.   Description                   Builder
    0x0010 : [0x0010, 200, "FSMD", "File system meta-data format", FSMD_200 ],
    0x0020 : [0x0020, 200, "RS+CFC8", "Reed-Solomon coded with Comma-free codewords", RS_CFC8_200 ],
    0x0100 : [0x0100, 200, "Dense", "Dense encoding", Dense_200 ],
