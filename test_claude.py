import anthropic

client = anthropic.Anthropic(api_key="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")

message = client.messages.create(
    model="claude-sonnet-4-6",
    max_tokens=1024,
    messages=[
        {"role": "user", "content": "Say hello and tell me one fact about HPLC validation."}
    ]
)

print(message.content[0].text)
