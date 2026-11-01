from openai import OpenAI

client = OpenAI()

response = client.responses.create(
    model="global.openai.gpt-5.6-sol",
    input="Can you explain the features of Amazon Bedrock?"
)
print(response)
