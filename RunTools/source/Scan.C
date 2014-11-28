void Scan(UInt_t event_id, const char* file_name)
{
    TFile f(file_name, "READ");
    std::ostringstream ss;
    ss << "EventId==" << event_id;
    const TSQLResult* result = events->Query("EventId", ss.str().c_str());
    if(result->GetRowCount() > 0)
        std::cerr << "found in " << file_name << std::endl;
}

